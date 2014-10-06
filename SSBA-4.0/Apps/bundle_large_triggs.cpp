#include "bundle_large_common.h"

#include "Base/v3d_timer.h"
#include "Math/v3d_nonlinlsq.h"

using namespace V3D;
using namespace std;

namespace
{

#define CAMERA_PARAM_TYPE 0
#define POINT_PARAM_TYPE 1

   struct SparseMetricBundleOptimizer;

   typedef InlineMatrix<double, 3, 6> Matrix3x6d;
   typedef InlineMatrix<double, 2, 6> Matrix2x6d;

   static int cameraParamDimensionFromMode(int mode)
   {
      switch (mode)
      {
         case FULL_BUNDLE_METRIC:       return 6;
         case FULL_BUNDLE_FOCAL_LENGTH: return 7;
         case FULL_BUNDLE_RADIAL:       return 9;
      }
      return 0;
   }

//**********************************************************************

#if !defined(USE_TUKEYS_BIWEIGHT)
   inline double rho(double const tau2, double const r2)
   {
      return (r2 < tau2) ? (r2*(1.0 - r2/2.0/tau2)/2.0) : (tau2/4.0);
   }

   inline double drho(double const tau2, double const r2)
   {
      return (r2 < tau2) ? (1.0 - r2/tau2)/2.0 : 0.0;
   }

   inline double ddrho(double const tau2, double const r2)
   {
      return (r2 < tau2) ? -1.0/tau2/2.0 : 0.0;
   }
#else
   inline double rho(double const tau2, double const r2)
   {
      double const r4 = r2*r2, tau4 = tau2*tau2;
      return (r2 < tau2) ? r2*(3.0 - 3*r2/tau2 + r4/tau4)/6.0f : tau2/6.0;
   }

   inline double drho(double const tau2, double const r2)
   {
      double const r4 = r2*r2, tau4 = tau2*tau2;
      return (r2 < tau2) ? (1.0 - 2*r2/tau2 + r4/tau4)/2.0 : 0.0;
   }

   inline double ddrho(double const tau2, double const r2)
   {
      double const tau4 = tau2*tau2;
      return (r2 < tau2) ? (r2/tau4 - 1.0/tau2) : 0.0;
   }
#endif

//**********************************************************************

   struct BundleCostFunction : public NLSQ_CostFunction
   {
         BundleCostFunction(int const mode, std::vector<int> const& usedParamTypes,
                            double const inlierThreshold,
                            vector<CameraMatrix> const& cams,
                            vector<SimpleDistortionFunction> const& distortions,
                            vector<Vector3d > const& Xs,
                            vector<Vector2d > const& measurements,
                            Matrix<int> const& correspondingParams)
            : NLSQ_CostFunction(usedParamTypes, correspondingParams, 2),
              _mode(mode), _cams(cams), _distortions(distortions), _Xs(Xs),  _measurements(measurements),
              _inlierThreshold(inlierThreshold), _sqrInlierThreshold(inlierThreshold*inlierThreshold)
         { }

         Vector2d projectPoint(Vector3d const& X, int i) const
         {
#if 0
            return _cams[i].projectPoint(_distortions[i], X);
#else
            Vector3d const XX = _cams[i].transformPointIntoCameraSpace(X);
            Vector2d const xu(XX[0] / XX[2],  XX[1] / XX[2]);
            Vector2d const xd = _distortions[i](xu);
            return _cams[i].getFocalLength() * xd;
#endif
         }

         virtual void evalResidual(int const k, Vector<double>& e) const
         {
            unsigned view  = _correspondingParams[k][0];
            unsigned point = _correspondingParams[k][1];

            Vector2d const q = this->projectPoint(_Xs[point], view);
            Vector2d const r = q - _measurements[k];

            e[0] = r[0];
            e[1] = r[1];
         }

         virtual double evalCost(VectorArray<double> const& residuals, Vector<double> const& W) const
         {
            double res = 0.0;

            for (int k = 0; k < _nMeasurements; ++k)
               res += rho(_sqrInlierThreshold, sqrNorm_L2(residuals[k]));

            return res;
         } // end evalCost()

         virtual double getWeight(Vector<double> const& r) const { return 1; }

         virtual void fillJacobian(int const whichParam, int const paramIx, int const k, Matrix<double>& Jdst) const
         {
            int const view  = _correspondingParams[k][0];
            int const point = _correspondingParams[k][1];

            Vector3d XX;
            Matrix3x6d dXX_dRT;
            Matrix3x3d dXX_dX;
            this->poseDerivatives(view, point, XX, dXX_dRT, dXX_dX);

            Vector2d xu; // undistorted image point
            xu[0] = XX[0] / XX[2];
            xu[1] = XX[1] / XX[2];

            Vector2d const xd = _distortions[view](xu); // distorted image point

            double const focalLength = _cams[view].getFocalLength();

            Matrix2x2d dp_dxd;
            dp_dxd[0][0] = focalLength; dp_dxd[0][1] = 0;
            dp_dxd[1][0] = 0;           dp_dxd[1][1] = focalLength;

            Matrix2x3d dxu_dXX;
            dxu_dXX[0][0] = 1.0f / XX[2]; dxu_dXX[0][1] = 0;            dxu_dXX[0][2] = -XX[0]/(XX[2]*XX[2]);
            dxu_dXX[1][0] = 0;            dxu_dXX[1][1] = 1.0f / XX[2]; dxu_dXX[1][2] = -XX[1]/(XX[2]*XX[2]);

            Matrix2x2d dxd_dxu = _distortions[view].derivativeWrtUndistortedPoint(xu);
            Matrix2x2d dp_dxu = dp_dxd * dxd_dxu;
            Matrix2x3d dp_dXX = dp_dxu * dxu_dXX;

            makeZeroMatrix(Jdst);

            if (whichParam == 0)
            {
               switch (_mode)
               {
                  case FULL_BUNDLE_RADIAL:
                  {
                     Matrix2x2d dxd_dk1k2 = _distortions[view].derivativeWrtRadialParameters(xu);
                     Matrix2x2d d_dk1k2 = dp_dxd * dxd_dk1k2;
                     copyMatrixSlice(d_dk1k2, 0, 0, 2, 2, Jdst, 0, 7);
                     // No break here!
                  }
                  case FULL_BUNDLE_FOCAL_LENGTH:
                  {
                     Jdst[0][6] = xd[0];
                     Jdst[1][6] = xd[1];
                  }
                  case FULL_BUNDLE_METRIC:
                  {
                     Matrix2x6d dp_dRT;
                     multiply_A_B(dp_dXX, dXX_dRT, dp_dRT);
                     copyMatrixSlice(dp_dRT, 0, 0, 2, 6, Jdst, 0, 0);
                  }
               } // end switch
            }
            else
            {
               // Jacobian w.r.t. 3d points
               multiply_A_B(dp_dXX, dXX_dX, Jdst);
            } // end if
         } // end fillJacobian()

         virtual void multiply_JtJ(double const lambda, int const k, Vector<double> const& residual, Matrix<double> const& J1, Matrix<double> const& J2,
                                   Matrix<double>& J1tJ2)
         {
            double const r2 = sqrNorm_L2(residual);
            double const rho_   = rho(_sqrInlierThreshold, r2);
            double const drho_  = drho(_sqrInlierThreshold, r2);
            double const ddrho_ = ddrho(_sqrInlierThreshold, r2);

            Matrix2x2d H;
            H[0][0] = H[1][1] = drho_;
            H[0][1] = H[1][0] = 0;

            if (drho_ + 2.0*ddrho_*r2 >= 0.0)
            {
               H[0][0] += 2.0*ddrho_*residual[0]*residual[0];
               H[1][1] += 2.0*ddrho_*residual[1]*residual[1];
               H[0][1] = H[1][0] = 2.0*ddrho_*residual[0]*residual[1];
            }
            // else
            // {
            //    double const s = drho_/std::max(1e-10, r2);
            //    H[0][0] -= s*residual[0]*residual[0];
            //    H[1][1] -= s*residual[1]*residual[1];
            //    H[0][1] = H[1][0] = -s*residual[0]*residual[1];
            // } // end if

            Matrix<double> H_J2(2, J2.num_cols());
            multiply_At_B(H, J2, H_J2);
            multiply_At_B(J1, H_J2, J1tJ2);
         }

         virtual void multiply_Jt_e(double const lambda, int const paramType, int const k, Matrix<double> const& Jk, Vector<double> const& residual,
                                    Vector<double>& Jkt_e)
         {
            double const r2 = sqrNorm_L2(residual);
            double const drho_ = drho(_sqrInlierThreshold, r2);

            Vector2d r;
            r[0] = drho_ * residual[0];
            r[1] = drho_ * residual[1];
            multiply_At_v(Jk, r, Jkt_e);
         }

      protected:
         void poseDerivatives(int i, int j, Vector3d& XX, Matrix3x6d& d_dRT, Matrix3x3d& d_dX) const
         {
            XX = _cams[i].transformPointIntoCameraSpace(_Xs[j]);

            // See Frank Dellaerts bundle adjustment tutorial.
            // d(dR * R0 * X + t)/d omega = -[R0 * X]_x
            Matrix3x3d J;
            makeCrossProductMatrix(XX - _cams[i].getTranslation(), J);
            scaleMatrixIP(-1.0, J);

            // Now the transformation from world coords into camera space is xx = Rx + T
            // Hence the derivative of x wrt. T is just the identity matrix.
            makeIdentityMatrix(d_dRT);
            copyMatrixSlice(J, 0, 0, 3, 3, d_dRT, 0, 3);

            // The derivative of Rx+T wrt x is just R.
            copyMatrix(_cams[i].getRotation(), d_dX);
         } // end poseDerivatives()

         int const _mode;

         vector<CameraMatrix>             const& _cams;
         vector<SimpleDistortionFunction> const& _distortions;
         vector<Vector3d>                 const& _Xs;

         double const _inlierThreshold, _sqrInlierThreshold;

         vector<Vector2d > const& _measurements;

         friend struct SparseMetricBundleOptimizer;
   }; // end struct BundleCostFunction

//**********************************************************************

   struct SparseMetricBundleOptimizer : public NLSQ_LM_Optimizer
   {
         typedef NLSQ_LM_Optimizer Base;

         SparseMetricBundleOptimizer(int const mode, NLSQ_ParamDesc const& paramDesc,
                                     std::vector<NLSQ_CostFunction *> const& costFunctions,
                                     double const inlierThreshold,
                                     vector<CameraMatrix>& cams,
                                     vector<SimpleDistortionFunction>& distortions,
                                     vector<Vector3d >& Xs)
            : Base(paramDesc, costFunctions), _mode(mode),
              _inlierThreshold(inlierThreshold), _sqrInlierThreshold(inlierThreshold*inlierThreshold),
              _cams(cams), _distortions(distortions), _Xs(Xs),
              _savedTranslations(cams.size()), _savedRotations(cams.size()), _savedFocalLengths(cams.size()), _savedDistortions(cams.size()),
              _savedXs(Xs.size()), _cachedParamLength(0.0)
         {
            // Since we assume that BA does not alter the inputs too much,
            // we compute the overall length of the parameter vector in advance
            // and return that value as the result of getParameterLength().
            for (int i = 0; i < _cams.size(); ++i)
            {
               _cachedParamLength += sqrNorm_L2(_cams[i].getTranslation());
               _cachedParamLength += 3.0; // Assume eye(3) for R.
            }

            if (mode >= FULL_BUNDLE_FOCAL_LENGTH)
               for (int i = 0; i < _cams.size(); ++i)
               {
                  double const f = _cams[i].getFocalLength();
                  _cachedParamLength += f*f;
               } // end for (i)

            for (int j = 0; j < _Xs.size(); ++j) _cachedParamLength += sqrNorm_L2(_Xs[j]);

            _cachedParamLength = sqrt(_cachedParamLength);
         }

         virtual double getParameterLength() const
         {
            return _cachedParamLength;
         }

         virtual void updateParameters(int const paramType, VectorArrayAdapter<double> const& delta)
         {
            switch (paramType)
            {
               case CAMERA_PARAM_TYPE:
               {
                  Vector3d T, omega;
                  Matrix3x3d R0, dR;

                  for (int i = 0; i < _cams.size(); ++i)
                  {
                     T = _cams[i].getTranslation();
                     T[0] += delta[i][0]; T[1] += delta[i][1]; T[2] += delta[i][2];
                     _cams[i].setTranslation(T);

                     // Create incremental rotation using Rodriguez formula.
                     R0 = _cams[i].getRotation();
                     omega[0] = delta[i][3]; omega[1] = delta[i][4]; omega[2] = delta[i][5];
                     createRotationMatrixRodrigues(omega, dR);
                     _cams[i].setRotation(dR * R0);

                     switch (_mode)
                     {
                        case FULL_BUNDLE_RADIAL:
                        {
                           _distortions[i].k1 += delta[i][7];
                           _distortions[i].k2 += delta[i][8];
                        }
                        case FULL_BUNDLE_FOCAL_LENGTH:
                        {
                           Matrix3x3d K = _cams[i].getIntrinsic();
                           K[0][0] += delta[i][6];
                           K[1][1] += delta[i][6];
                           _cams[i].setIntrinsic(K);
                        }
                     } // end switch (_mode)
                  }
                  break;
               }
               case POINT_PARAM_TYPE:
               {
                  for (int j = 0; j < _Xs.size(); ++j)
                  {
                     _Xs[j][0] += delta[j][0];
                     _Xs[j][1] += delta[j][1];
                     _Xs[j][2] += delta[j][2];
                  }
                  break;
               }
               default:
                  assert(false);
            } // end switch (paramType)
         } // end updateParametersA()

         virtual void saveAllParameters()
         {
            for (int i = 0; i < _cams.size(); ++i)
            {
               _savedTranslations[i] = _cams[i].getTranslation();
               _savedRotations[i]    = _cams[i].getRotation();
               _savedFocalLengths[i] = _cams[i].getFocalLength();
            }
            _savedXs = _Xs;
            _savedDistortions = _distortions;
         }

         virtual void restoreAllParameters()
         {
            for (int i = 0; i < _cams.size(); ++i)
            {
               _cams[i].setTranslation(_savedTranslations[i]);
               _cams[i].setRotation(_savedRotations[i]);
               Matrix3x3d K = _cams[i].getIntrinsic();
               K[0][0] = K[1][1] = _savedFocalLengths[i];
               _cams[i].setIntrinsic(K);
            }
            _Xs = _savedXs;
            _distortions = _savedDistortions;
         }

      protected:
         int const _mode;
         double const _inlierThreshold, _sqrInlierThreshold;

         vector<CameraMatrix>& _cams;
         vector<Vector3d >&    _Xs;
         vector<SimpleDistortionFunction>& _distortions;

         vector<Vector3d >   _savedTranslations;
         vector<Matrix3x3d > _savedRotations;
         vector<Vector3d >   _savedXs;

         vector<double>                   _savedFocalLengths;
         vector<SimpleDistortionFunction> _savedDistortions;

         double _cachedParamLength;
   }; // end struct SparseMetricBundleOptimizer

//**********************************************************************

   int
   adjustStructureAndMotion(int const mode, 
                            vector<CameraMatrix>& cams,
                            vector<SimpleDistortionFunction>& distortions,
                            vector<Vector3d >& Xs,
                            vector<Vector2d > const& measurements2d,
                            vector<int> const& correspondingView,
                            vector<int> const& correspondingPoint,
                            double inlierThreshold)
   {
      NLSQ_ParamDesc paramDesc;
      paramDesc.nParamTypes = 2;
      paramDesc.dimension[CAMERA_PARAM_TYPE] = cameraParamDimensionFromMode(mode);
      paramDesc.dimension[POINT_PARAM_TYPE]  = 3;

      paramDesc.count[CAMERA_PARAM_TYPE] = cams.size();
      paramDesc.count[POINT_PARAM_TYPE]  = Xs.size();

      vector<int> usedParamTypes;
      usedParamTypes.push_back(CAMERA_PARAM_TYPE);
      usedParamTypes.push_back(POINT_PARAM_TYPE);

      Matrix<int> correspondingParams(measurements2d.size(), paramDesc.nParamTypes);
      for (int k = 0; k < correspondingParams.num_rows(); ++k)
      {
         correspondingParams[k][0] = correspondingView[k];
         correspondingParams[k][1] = correspondingPoint[k];
      }

      BundleCostFunction costFun(mode, usedParamTypes, inlierThreshold, cams, distortions, Xs, measurements2d, correspondingParams);
      vector<NLSQ_CostFunction *> costFunctions;
      costFunctions.push_back(&costFun);

      SparseMetricBundleOptimizer opt(mode, paramDesc, costFunctions, inlierThreshold, cams, distortions, Xs);
      opt.updateThreshold = 1e-12;
      //opt.updateThreshold = 1e-20;
      opt.maxIterations = 100; //params.nIterations;
      opt.tau = 1e-3;

      Timer t("BA");
      t.start();
      opt.minimize();
      t.stop();
      cout << "Time per iteration: " << t.getTime() / opt.currentIteration << endl;

      //params.lambda = opt.lambda;
      return opt.status;
   }

} // end namespace <>

int
main(int argc, char * argv[])
{
   if (argc != 2)
   {
      cerr << "Usage: " << argv[0] << " <sparse reconstruction file>" << endl;
      return -1;
   }

   ifstream is(argv[1]);
   if (!is)
   {
      cerr << "Cannot open " << argv[1] << endl;
      return -2;
   }

   double const avg_focal_length = AVG_FOCAL_LENGTH;
   V3D::optimizerVerbosenessLevel = 1;

   cout.precision(10);

   int N, M, K;
   is >> N >> M >> K;
   cout << "N (cams) = " << N << " M (points) = " << M << " K (measurements) = " << K << endl;

   cout << "Reading image measurements..." << endl;
   vector<Vector2d> measurements(K);
   vector<int> correspondingView(K, -1);
   vector<int> correspondingPoint(K, -1);
   for (int k = 0; k < K; ++k)
   {
      is >> correspondingView[k];
      is >> correspondingPoint[k];
      is >> measurements[k][0] >> measurements[k][1];
      measurements[k][0] /= avg_focal_length;
      measurements[k][1] /= avg_focal_length;
   } // end for (k)
   cout << "Done." << endl;

   cout << "Reading cameras..." << endl;
   vector<CameraMatrix> cams(N);
   vector<SimpleDistortionFunction> distortions(N);
   for (int i = 0; i < N; ++i)
   {
      Vector3d om, T;
      double f, k1, k2;
      is >> om[0] >> om[1] >> om[2];
      is >> T[0] >> T[1] >> T[2];
      is >> f >> k1 >> k2;

      Matrix3x3d K; makeIdentityMatrix(K);
      K[0][0] = K[1][1] = -f / avg_focal_length;
      cams[i].setIntrinsic(K);
      cams[i].setTranslation(T);

      Matrix3x3d R;
      createRotationMatrixRodrigues(om, R);
      cams[i].setRotation(R);

      double const f2 = f*f;
      distortions[i].k1 = k1 * f2;
      distortions[i].k2 = k2 * f2 * f2;

      //cout << "k1 = " << k1 << " k2 = " << k2 << endl;
   } // end for (i)
   cout << "Done." << endl;

   cout << "Reading 3D point..." << endl;
   vector<Vector3d > Xs(M);
   for (int j = 0; j < M; ++j) is >> Xs[j][0] >> Xs[j][1] >> Xs[j][2];
   cout << "Done." << endl;

   double init_ratio = showErrorStatistics(avg_focal_length, inlier_threshold, cams, distortions, Xs, measurements, correspondingView, correspondingPoint);
   double const E_init = showObjective(avg_focal_length, inlier_threshold, cams, distortions, Xs, measurements, correspondingView, correspondingPoint);

   for (int i = 0; i < 10; ++i) cout << "f[" << i << "] = " << cams[i].getFocalLength() << endl;

   adjustStructureAndMotion(bundle_mode, cams, distortions, Xs, measurements, correspondingView, correspondingPoint,
                            inlier_threshold/avg_focal_length);

   for (int i = 0; i < 10; ++i) cout << "f[" << i << "] = " << cams[i].getFocalLength() << endl;

   double final_ratio = showErrorStatistics(avg_focal_length, inlier_threshold, cams, distortions, Xs, measurements, correspondingView, correspondingPoint);
   //showErrorStatistics(KMat, cams, Xs, measurements, correspondingView, correspondingPoint);
   double const E_final = showObjective(avg_focal_length, inlier_threshold, cams, distortions, Xs, measurements, correspondingView, correspondingPoint);
   cout << "E_init = " << E_init << " E_final = " << E_final << " initial ratio = " << init_ratio << " final ratio = " << final_ratio << endl;

   return 0;
}
