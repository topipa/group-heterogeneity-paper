#include <Rcpp.h>
#include <stan/math.hpp>
#include "stan/math/mix/mat.hpp"
#include <RcppEigen.h>

using namespace Rcpp;

class sqr_exponential_ard_kernel_cpp;
RCPP_EXPOSED_CLASS(sqr_exponential_ard_kernel_cpp);

// [[Rcpp::plugins(cpp14)]]

template<typename T>
struct sqr_exp_ard_kernel_alt_functor {
  const Eigen::Matrix<double, Eigen::Dynamic, 1> phi_;
  const Eigen::Matrix<double, Eigen::Dynamic, 1> phi2_;
  Eigen::Matrix<T, Eigen::Dynamic, 1> x_i_;
  Eigen::Matrix<T, Eigen::Dynamic, 1> x_i_scaled_;
  T x_i_squared_;
  sqr_exp_ard_kernel_alt_functor(const Eigen::Matrix<double, Eigen::Dynamic, 1>& phi,
                                 const Eigen::Matrix<double, Eigen::Dynamic, 1>& phi2,
                             Eigen::Matrix<T, Eigen::Dynamic, 1> x_i) : phi_(phi), phi2_(phi2), x_i_(x_i) {

    Eigen::Matrix<double, Eigen::Dynamic, 1> lengthscale = phi2_;
    lengthscale = lengthscale.cwiseInverse();

    x_i_scaled_ = x_i_.cwiseProduct(lengthscale);

    x_i_squared_ = (x_i_scaled_.transpose() * x_i_scaled_).sum();
  }
  void update_x_i(Eigen::Matrix<T, Eigen::Dynamic, 1> x_i) {
    Eigen::Matrix<double, Eigen::Dynamic, 1> lengthscale = phi2_;
    lengthscale = lengthscale.cwiseInverse();
    x_i_ = x_i;
    x_i_scaled_ = x_i_.cwiseProduct(lengthscale);
    x_i_squared_ = (x_i_scaled_.transpose() * x_i_scaled_).sum();
  }
  template<typename G>
  G operator() (Eigen::Matrix<G, Eigen::Dynamic, 1> x_j) const {
    double alpha = phi_(0);
    Eigen::Matrix<double, Eigen::Dynamic, 1> lengthscale = phi2_;
    lengthscale = lengthscale.cwiseInverse();

    alpha = alpha * alpha;

    Eigen::Matrix<G, Eigen::Dynamic, 1> x_j_scaled;

    x_j_scaled = x_j.cwiseProduct(lengthscale);

    G x_j_squared;

    x_j_squared = (x_j_scaled.transpose() * x_j_scaled).sum();

    G k_ij = (x_i_squared_ + x_j_squared
                - 2 * (x_i_scaled_.transpose() * x_j_scaled).sum());
    k_ij = alpha * stan::math::exp(-0.5 * k_ij);

    return k_ij;
  }
};


template <typename T>
struct sqr_exp_ard_kernel_full_functor {
  const Eigen::Matrix<double, Eigen::Dynamic, 1> phi_;
  const Eigen::Matrix<double, Eigen::Dynamic, 1> phi2_;
  const int N_;
  const int M_;
  const int D_;
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> X_;
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Z_;
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> K_;
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> G_;
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> H_;

  sqr_exp_ard_kernel_full_functor(const Eigen::Matrix<double, Eigen::Dynamic, 1>& phi,
                                  const Eigen::Matrix<double, Eigen::Dynamic, 1>& phi2,
                              Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& X,
                              Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& Z) :
    phi_(phi), phi2_(phi2), N_(X.rows()), M_(Z.rows()), D_(X.cols()), X_(X), Z_(Z), \
    H_(Eigen::MatrixXd::Zero(N_*M_, D_*D_)), G_(Eigen::MatrixXd::Zero(N_*M_, D_)), \
    K_(Eigen::MatrixXd::Zero(N_, M_))  { }

  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> operator() () {
    Eigen::Matrix<T, 1, Eigen::Dynamic> x_i = X_.row(0);
    sqr_exp_ard_kernel_alt_functor<T> k_ij(phi_, phi2_, x_i.transpose());

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> K(N_, M_);
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> G(N_*M_, D_);
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> H(N_*M_, D_*D_);
    K = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(N_, M_);
    for (int i = 0; i < N_; ++i) {
      if (i > 0) {
        Eigen::Matrix<T, 1, Eigen::Dynamic> x_i = X_.row(i);
        k_ij.update_x_i(x_i.transpose());
      }
      for (int j = 0; j < M_; ++j) {
        Eigen::Matrix<T, Eigen::Dynamic, 1> z_j = Z_.row(j).transpose();
        Eigen::Matrix<T, Eigen::Dynamic, 1> grad;
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> h;
        double f;
        stan::math::hessian(k_ij, z_j, f, grad, h);
        K(i,j) = f;

        for (int d = 0; d < D_; ++d) {
          G(i * M_ + j, d) = grad(d);
          for (int k = 0; k < D_; ++k)
            H(i * M_ + j, d * D_ + k) = h(d, k);
        }
      }
    }

    G_ = G;
    H_ = H;

    return K;
  }
};

class sqr_exponential_ard_kernel_cpp {
public:
  sqr_exponential_ard_kernel_cpp(Eigen::Matrix<double, Eigen::Dynamic, 1> phi,
                                 Eigen::Matrix<double, Eigen::Dynamic, 1> phi2,
                             Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> X,
                             Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Z) {
    sqr_exp_ard_kernel_full_functor<double> K(phi,phi2, X, Z);

    Kernel = K();
    dK_dX = K.G_;
    d2K_dXdX = K.H_;
  }
  ~sqr_exponential_ard_kernel_cpp() { }

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Kernel;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dK_dX;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> d2K_dXdX;
};

RCPP_MODULE(sqr_exponential_ard_kernel){
  class_<sqr_exponential_ard_kernel_cpp>( "sqr_exponential_ard_kernel" )

  .constructor<Eigen::Matrix<double, Eigen::Dynamic, 1>,
  Eigen::Matrix<double, Eigen::Dynamic, 1>,
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>,
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>()

  .field("Kernel", &sqr_exponential_ard_kernel_cpp::Kernel)
  .field("dK_dX", &sqr_exponential_ard_kernel_cpp::dK_dX)
  .field("d2K_dXdX", &sqr_exponential_ard_kernel_cpp::d2K_dXdX)
  ;
}
