#include <iostream>
#include <cmath>

#include "vice.h"

template <typename T>
class TrapFunctor {
  private:
    T depth_, duration_, ingress_;
    long int n_eval_;

  public:

    TrapFunctor (T depth, T duration, T ingress)
      : depth_(depth), duration_(duration), ingress_(ingress), n_eval_(0) {}

    void reset () { n_eval_ = 0; }
    long int get_n_eval () const { return n_eval_; }
    T get_full_duration () const { return duration_ + ingress_; }

    T operator() (T t) {
      n_eval_++;
      t = std::abs(t);
      if (t < duration_) {
        return 1 - depth_;
      } else if (t < duration_ + ingress_) {
        return 1 + (t - duration_ - ingress_) * depth_ / ingress_;
      }
      return 1;
    }

    T integrate (T lower, T upper) {
      T dt = duration_ + ingress_;
      T integral = T(0);

      if (lower < -duration_-ingress_) {
        T b = std::min(-duration_-ingress_, upper);
        integral += b - lower;
        lower = b;
      }

      if (lower < -duration_) {
        T b = std::min(-duration_, upper);

        integral += (-b*b*depth_/2 - b*(depth_*duration_ + depth_*ingress_ - ingress_) + depth_*lower*lower/2 + lower*(depth_*duration_ + depth_*ingress_ - ingress_))/ingress_;

        lower = b;
      }

      if (lower < duration_) {
        T b = std::min(duration_, upper);
        integral += (1 - depth_) * (b - lower);
        lower = b;
      }

      if (lower < duration_ + ingress_) {
        T b = std::min(duration_ + ingress_, upper);

        integral += (b*b*depth_/2 - b*(depth_*duration_ + depth_*ingress_ - ingress_) - depth_*lower*lower/2 + lower*(depth_*duration_ + depth_*ingress_ - ingress_))/ingress_;

        lower = b;
      }

      if (lower < upper) {
        integral += upper - lower;
      }

      return integral;
    }
};


template <typename Scalar, typename Functor, typename Integrator, typename ParamType>
std::tuple<std::vector<Scalar>, long int> benchmark_integrate (Functor& func, Integrator& integrator, const std::vector<Scalar>& time, const Scalar texp, ParamType param) {
  func.reset();
  std::vector<Scalar> flux(time.size());

  for (size_t i = 0; i < time.size(); ++i) {
    Scalar t = time[i];
    Scalar lower = t - 0.5*texp;
    Scalar upper = t + 0.5*texp;
    flux[i] = integrator(func, lower, upper, param) / texp;
  }
  return std::make_tuple(flux, func.get_n_eval());
}


int main () {
  double texp = 0.1;

  double depth    = 0.1;
  double duration = 0.1;
  double ingress  = 0.2;
  TrapFunctor<double> func(depth, duration, ingress);

  double tmin = -duration - 2*ingress;
  double tmax =  duration + 2*ingress;
  int npts = 500;
  double dt = (tmax - tmin) / (npts - 1);

  int neval = 15;
  double tol = 1e-7;

  std::cout << "time,flux,truth,reimann,n_reimann,trap_fixed,n_trap_fixed,simp_fixed,n_simp_fixed,trap_adapt,n_trap_adapt,simp_adapt,n_simp_adapt,\n";
  std::cout.precision(25);
  std::cout << std::scientific;

  for (int i = 0; i < npts; ++i) {
    double t = tmin + i * dt;
    double lower = t - 0.5*texp;
    double upper = t + 0.5*texp;

    func.reset();
    double flux = func(t);
    double truth = func.integrate(lower, upper) / texp;

    std::cout << t << "," << flux << "," << truth << ",";

    func.reset();
    double f_reimann = confessional::reimann::integrate_fixed(func, lower, upper, neval) / texp;
    std::cout << f_reimann << "," << func.get_n_eval() << ",";

    func.reset();
    double f_trap_fixed = confessional::trapezoid::integrate_fixed(func, lower, upper, neval) / texp;
    std::cout << f_trap_fixed << "," << func.get_n_eval() << ",";

    func.reset();
    double f_simp_fixed = confessional::simpson::integrate_fixed(func, lower, upper, neval) / texp;
    std::cout << f_simp_fixed << "," << func.get_n_eval() << ",";

    func.reset();
    double f_trap_adapt = confessional::trapezoid::integrate_adapt(func, lower, upper, tol, 50) / texp;
    std::cout << f_trap_adapt << "," << func.get_n_eval() << ",";

    func.reset();
    double f_simp_adapt = confessional::simpson::integrate_adapt(func, lower, upper, tol, 50) / texp;
    std::cout << f_simp_adapt << "," << func.get_n_eval() << ",";

    std::cout << "\n";
  }
}
