#include <iostream>
#include <cmath>

#include <Eigen/Core>

#include "vice/integrate.h"
#include "vice/functors.h"


int main () {
  double texp = 0.1;

  //double depth    = 0.1;
  //double duration = 0.1;
  //double ingress  = 0.2;
  //vice::functors::TrapFunctor<double> func(depth, duration, ingress);

  double tau = 0.5;
  double r = 0.1;
  double b = 0.5;
  Eigen::VectorXd u(2);
  u << 0.3, 0.2;
  vice::functors::StarryFunctor<double> func(u, r, b, tau);

  //double tmin = -duration - 2*ingress;
  //double tmax =  duration + 2*ingress;
  double tmin = -tau - texp;
  double tmax =  tau + texp;
  int npts = 500;
  double dt = (tmax - tmin) / (npts - 1);

  int neval = 15;
  double tol = 1e-9;

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
    double f_reimann = vice::reimann::integrate_fixed(func, lower, upper, neval) / texp;
    std::cout << f_reimann << "," << func.get_n_eval() << ",";

    func.reset();
    double f_trap_fixed = vice::trapezoid::integrate_fixed(func, lower, upper, neval) / texp;
    std::cout << f_trap_fixed << "," << func.get_n_eval() << ",";

    func.reset();
    double f_simp_fixed = vice::simpson::integrate_fixed(func, lower, upper, neval) / texp;
    std::cout << f_simp_fixed << "," << func.get_n_eval() << ",";

    func.reset();
    double f_trap_adapt = vice::trapezoid::integrate_adapt(func, lower, upper, tol, 50) / texp;
    std::cout << f_trap_adapt << "," << func.get_n_eval() << ",";

    func.reset();
    double f_simp_adapt = vice::simpson::integrate_adapt(func, lower, upper, tol, 50) / texp;
    std::cout << f_simp_adapt << "," << func.get_n_eval() << ",";

    std::cout << "\n";
  }
}
