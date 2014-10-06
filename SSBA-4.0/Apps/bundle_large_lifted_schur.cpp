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
   inline double kappa(double const tau, double const w2) { return 0.70710678118*tau*(w2 - 1); }
   inline double dkappa(double const tau, double const w2) { return 0.70710678118*tau; }
#else
   inline double kappa(double const tau, double const w2)
   {
      double const f = sqrt(1.0/3);
      double const w = sqrt(w2);
      return f * tau * (w - 1)*sqrt(2*w+1);
   }
   inline double dkappa(double const tau, double const w2)
   {
      double const f = sqrt(3.0)/2;
      double const w = sqrt(w2);
      return f * tau / sqrt(2*w+1);
   }
#endif

// The mapping w |-> h(w)
#if 1
   inline double omega(double const w) { return w; }
   inline double domega(double const w) { return 1.0; }
   inline double omega2_inv(double const w) { return sqrt(w); }
#elif 0
   //double const slope_w = 0.25; //3.0;
   double const slope_w = 4.0;
   inline double omega(double const w) { return slope_w*w; }
   inline double domega(double const w) { return slope_w; }
   inline double omega2_inv(double const w) { return sqrt(w)/slope_w; }
#elif 1
   inline double omega(double const w) { return 0.5*(w/sqrt(1+w*w)+1); }
   inline double domega(double const w) { return 0.5*sqrt(1+w*w)/(1+w*w); }
   inline double omega2_inv(double const w) { double const ww = sqrt(w); return (2*ww-1.0)/sqrt(1.0 - sqr(2*ww-1)); }
// #elif 1
//    inline double omega(double const w) { return w*w; }
//    inline double domega(double const w) { return 2.0*w; }
//    inline double omega2_inv(double const w) { return sqrt(sqrt(w)); }
// #elif 0
//    inline double omega(double const w) { return sqr(std::max(0.0, w)); }
//    inline double domega(double const w) { return (w > 0) ? 2.0*w : 0.0; }
//    inline double omega2_inv(double const w) { return sqrt(sqrt(w)); }
#elif 1
   // h(w) = exp(w)
   double const exp_const = 1.0/4.0;
   inline double omega(double const w) { return exp(w * exp_const); }
   inline double domega(double const w) { return exp_const * exp(w * exp_const); }
   inline double omega2_inv(double const w) { return log(w)/2.0/exp_const; }
#else
   // h(w) = soft-max(0, w), logistic
   double const exp_const = 4.0;
   inline double omega(double const w) { return exp_const * log(1.0 + exp(w/exp_const)); }
   inline double domega(double const w) { double const exp_w = exp(w/exp_const); return exp_w/(1.0 + exp_w); }
   inline double omega2_inv(double const w) { return exp_const * log(exp(sqrt(w)/exp_const) - 1.0); }
#endif


//**********************************************************************

   struct BundleCostFunction : public NLSQ_CostFunction
   {
         BundleCostFunction(int const mode, std::vector<int> const& usedParamTypes,
                            double const inlierThreshold,
                            vector<CameraMatrix> const& cams,
                            vector<SimpleDistortionFunction> const& distortions,
                            vector<Vector3d > const& Xs,
                            vector<double> const& weights,
                            vector<Vector2d > const& measurements,
                            Matrix<int> const& correspondingParams)
            : NLSQ_CostFunction(usedParamTypes, correspondingParams, 2),
              _mode(mode), _cams(cams), _distortions(distortions), _Xs(Xs), _weights(weights), _measurements(measurements),
              _inlierThreshold(inlierThreshold), _sqrInlierThreshold(inlierThreshold*inlierThreshold),
              _Jt_e_cams(measurements.size(), cameraParamDimensionFromMode(mode)), _Jt_e_Xs(measurements.size(), 3)
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
#if 0
            double res = 0.0;

            for (int k = 0; k < _nMeasurements; ++k)
            {
               double const w2 = sqr(omega(_weights[k]));
               res += psi_hat(_sqrInlierThreshold, sqrNorm_L2(residuals[k]), w2);
            }

            return res/2;
#else
            double E_data = 0.0, E_reg = 0.0;

            for (int k = 0; k < _nMeasurements; ++k)
            {
               double const w2 = sqr(omega(_weights[k]));
               E_data += w2*sqrNorm_L2(residuals[k]);
               E_reg += sqr(kappa(_inlierThreshold, w2));
            }
            E_data /= 2.0; E_reg /= 2.0;
            //cout << "evalCost(): E_data = " << E_data << " E_reg = " << E_reg << endl;

            return E_data + E_reg;
#endif
         } // end evalCost()

         virtual double getWeight(Vector<double> const& r) const { return 1; }

         virtual void fillJacobian(int const whichParam, int const paramIx, int const k, Matrix<double>& Jdst) const
         {
            int const view  = _correspondingParams[k][0];
            int const point = _correspondingParams[k][1];

            //double const w = _weights[k];

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

            //scaleMatrixIP(_weights[k], Jdst); // Actually used Jacobian is W*J
         } // end fillJacobian()

         virtual void multiply_JtJ(double const lambda, int const k, Vector<double> const& residual, Matrix<double> const& J1, Matrix<double> const& J2,
                                   Matrix<double>& J1tJ2)
         {
            double const w2 = sqr(omega(_weights[k])), dw2 = sqr(domega(_weights[k]));

            Matrix2x2d C;
            makeOuterProductMatrix(residual, C);
#if 1
            scaleMatrixIP(-1.0/(sqrNorm_L2(residual) + 4.0*w2*sqr(dkappa(_inlierThreshold, w2)) + lambda/dw2), C);
#else
            scaleMatrixIP(-1.0/(sqrNorm_L2(residual) + 2.0*w2*_sqrInlierThreshold + lambda), C);
#endif
            C[0][0] += 1; C[1][1] += 1;
            scaleMatrixIP(w2, C);

            Matrix<double> C_J2(2, J2.num_cols());
            multiply_At_B(C, J2, C_J2);
            multiply_At_B(J1, C_J2, J1tJ2);
         }

         virtual void multiply_Jt_e(double const lambda, int const paramType, int const k, Matrix<double> const& Jk, Vector<double> const& residual,
                                    Vector<double>& Jkt_e)
         {
            double const r2 = sqrNorm_L2(residual), w2 = sqr(omega(_weights[k])), dw2 = sqr(domega(_weights[k]));
#if 1
            double const num = r2 + 2.0*dkappa(_inlierThreshold, w2)*kappa(_inlierThreshold, w2);
            double const denom = r2 + 4.0*w2*sqr(dkappa(_inlierThreshold, w2)) + lambda/dw2;
            double const f = w2 * (1.0 - num/denom);
#else
            double const num = lambda + _sqrInlierThreshold*(w2 + 1.0);
            double const denom = r2 + 2.0*_sqrInlierThreshold*w2 + lambda;
            double const f = w2 * num/denom;
#endif

            Vector2d r;
            r[0] = f * residual[0];
            r[1] = f * residual[1];
            multiply_At_v(Jk, r, Jkt_e);

            if (paramType == 0)
               multiply_At_v(Jk, residual, _Jt_e_cams[k]);
            else
               multiply_At_v(Jk, residual, _Jt_e_Xs[k]);
         }

         virtual bool forbid_derivative_caching() const { return true; }

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
         vector<double>                   const& _weights;

         double const _inlierThreshold, _sqrInlierThreshold;

         vector<Vector2d > const& _measurements;

         VectorArray<double> _Jt_e_cams, _Jt_e_Xs;

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
                                     vector<Vector3d >& Xs,
                                     vector<double>& weights)
            : Base(paramDesc, costFunctions), _mode(mode),
              _inlierThreshold(inlierThreshold), _sqrInlierThreshold(inlierThreshold*inlierThreshold),
              _cams(cams), _distortions(distortions), _Xs(Xs), _weights(weights),
              _savedTranslations(cams.size()), _savedRotations(cams.size()), _savedFocalLengths(cams.size()), _savedDistortions(cams.size()),
              _savedXs(Xs.size()), _cachedParamLength(0.0),
              _delta_cams(cams.size(), paramDesc.dimension[CAMERA_PARAM_TYPE]),
              _delta_Xs(Xs.size(), paramDesc.dimension[POINT_PARAM_TYPE])
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
                     copyVector(delta[i], _delta_cams[i]);

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
                     copyVector(delta[j], _delta_Xs[j]);

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

         virtual void finishUpdateParameters()
         {
            // Update the weights according to the Schur complement
            NLSQ_Residuals const& residuals = *_residuals[0];
            BundleCostFunction const& costFun = *(BundleCostFunction *)_costFunctions[0];

            Matrix<int> const& correspondingParams = costFun.correspondingParams();

            for (int k = 0; k < _weights.size(); ++k)
            {
               unsigned const i = correspondingParams[k][0];
               unsigned const j = correspondingParams[k][1];

               double const w = omega(_weights[k]), w2 = w*w, dw = domega(_weights[k]), dw2 = dw*dw;
               double const r2 = sqrNorm_L2(residuals._residuals[k]);

               double const rt_J_delta_cam = innerProduct(costFun._Jt_e_cams[k], _delta_cams[i]);
               double const rt_J_delta_Xs  = innerProduct(costFun._Jt_e_Xs[k], _delta_Xs[j]);
#if 1
               double const num = 2.0*dkappa(_inlierThreshold, w2)*kappa(_inlierThreshold, w2) + r2 + rt_J_delta_cam + rt_J_delta_Xs;
               double const denom = 4.0*w2*sqr(dkappa(_inlierThreshold, w2)) + r2 + this->lambda/dw2;
#else
               double const num = _sqrInlierThreshold*(w2-1) + r2 + rt_J_delta_cam + rt_J_delta_Xs;
               double const denom = 2.0*_sqrInlierThreshold*w2 + r2 + this->lambda;
#endif
               double const delta_w = -w/dw*num/denom;
               _weights[k] += delta_w;
            } // end for (k)
         } // end finishUpdateParameters()

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
            _savedWeights = _weights;
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
            _weights = _savedWeights;
         }

      protected:
         int const _mode;
         double const _inlierThreshold, _sqrInlierThreshold;

         vector<CameraMatrix>& _cams;
         vector<Vector3d >&    _Xs;
         vector<SimpleDistortionFunction>& _distortions;
         vector<double>&       _weights;

         vector<Vector3d >   _savedTranslations;
         vector<Matrix3x3d > _savedRotations;
         vector<Vector3d >   _savedXs;
         vector<double>      _savedWeights;

         vector<double>                   _savedFocalLengths;
         vector<SimpleDistortionFunction> _savedDistortions;

         double _cachedParamLength;

         VectorArray<double> _delta_cams, _delta_Xs;
   }; // end struct SparseMetricBundleOptimizer

//**********************************************************************

   int
   adjustStructureAndMotion(int const mode, 
                            vector<CameraMatrix>& cams,
                            vector<SimpleDistortionFunction>& distortions,
                            vector<Vector3d >& Xs,
                            vector<double>& weights,
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

      BundleCostFunction costFun(mode, usedParamTypes, inlierThreshold, cams, distortions, Xs, weights, measurements2d, correspondingParams);
      vector<NLSQ_CostFunction *> costFunctions;
      costFunctions.push_back(&costFun);

      SparseMetricBundleOptimizer opt(mode, paramDesc, costFunctions, inlierThreshold, cams, distortions, Xs, weights);
      opt.updateThreshold = 1e-12;
      //opt.updateThreshold = 1e-20;
      opt.maxIterations = 100; //params.nIterations;
      opt.tau = 1e-3;

      Timer t("BA");
      t.start();
      opt.minimize();
      t.stop();
      cout << "Time per iteration: " << t.getTime() / opt.currentIteration << endl;

      if (0)
      {
         // Reset the weights and restart NSLQ
         cout << "2nd round of NSQL" << endl;

         double const tau2 = sqr(inlierThreshold);

         for (int k = 0; k < weights.size(); ++k)
         {
#if 0
            weights[k] = omega2_inv(0.5);
#else
            int const i = correspondingView[k];
            int const j = correspondingPoint[k];
            Vector2d p = cams[i].projectPoint(distortions[i], Xs[j]);

            double const r2 = sqrNorm_L2(p - measurements2d[k]);
            weights[k] = omega2_inv(psi_weight(tau2, r2));
#endif
         } // end for (k)
         opt.tau = 1e-3;
         opt.minimize();
      }

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

   vector<double> weights(measurements.size(), omega2_inv(1.0));

   adjustStructureAndMotion(bundle_mode, cams, distortions, Xs, weights, measurements, correspondingView, correspondingPoint,
                            inlier_threshold/avg_focal_length);

   for (int i = 0; i < 10; ++i) cout << "f[" << i << "] = " << cams[i].getFocalLength() << endl;

   double final_ratio = showErrorStatistics(avg_focal_length, inlier_threshold, cams, distortions, Xs, measurements, correspondingView, correspondingPoint);
   //showErrorStatistics(KMat, cams, Xs, measurements, correspondingView, correspondingPoint);
   double const E_final = showObjective(avg_focal_length, inlier_threshold, cams, distortions, Xs, measurements, correspondingView, correspondingPoint);
   cout << "E_init = " << E_init << " E_final = " << E_final << " initial ratio = " << init_ratio << " final ratio = " << final_ratio << endl;

   return 0;
}
