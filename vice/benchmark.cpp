#include <tuple>

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include <Eigen/Core>

#include "vice/functors.h"
#include "vice/integrate.h"
#include "vice/benchmark.h"

namespace py = pybind11;

//
// Trapezoid functions
//
template <typename Scalar, typename ParamType, typename Integrator>
std::tuple<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>, Eigen::Matrix<long int, Eigen::Dynamic, 1>, double>
benchmark_trapezoid (Scalar depth, Scalar duration, Scalar ingress,
                     Eigen::Ref<const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>> time, const Scalar texp,
                     ParamType param)
{
  vice::functors::TrapFunctor<Scalar> func(depth, duration, ingress);
  auto integrator = Integrator();
  return vice::benchmark::integrate(func, integrator, time, texp, param);
}

template <typename Scalar>
Eigen::Matrix<Scalar, Eigen::Dynamic, 1>
trapezoid_exact (Scalar depth, Scalar duration, Scalar ingress,
                 Eigen::Ref<const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>> time, Scalar texp)
{
  Eigen::Matrix<Scalar, Eigen::Dynamic, 1> fluence(time.rows());
  vice::functors::TrapFunctor<Scalar> func(depth, duration, ingress);
  for (int i = 0; i < time.rows(); ++i) {
    Scalar lower = time(i) - 0.5*texp;
    Scalar upper = time(i) + 0.5*texp;
    fluence(i) = func.integrate(lower, upper);
  }
  return fluence;
}

template <typename Scalar>
Eigen::Matrix<Scalar, Eigen::Dynamic, 1>
trapezoid_flux (Scalar depth, Scalar duration, Scalar ingress,
                Eigen::Ref<const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>> time)
{
  Eigen::Matrix<Scalar, Eigen::Dynamic, 1> flux(time.rows());
  vice::functors::TrapFunctor<Scalar> func(depth, duration, ingress);
  for (int i = 0; i < time.rows(); ++i) {
    flux(i) = func(time(i));
  }
  return flux;
}


//
// Transit functions
//
template <typename Scalar, typename ParamType, typename Integrator>
std::tuple<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>, Eigen::Matrix<long int, Eigen::Dynamic, 1>, double>
benchmark_transit (Eigen::Ref<const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>> u, Scalar r, Scalar b, Scalar tau,
                   Eigen::Ref<const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>> time, const Scalar texp,
                   ParamType param)
{
  vice::functors::StarryFunctor<Scalar> func(u, r, b, tau);
  auto integrator = Integrator();
  return vice::benchmark::integrate(func, integrator, time, texp, param);
}

template <typename Scalar>
Eigen::Matrix<Scalar, Eigen::Dynamic, 1>
transit_exact (Eigen::Ref<const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>> u, Scalar r, Scalar b, Scalar tau,
               Eigen::Ref<const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>> time, Scalar texp)
{
  Eigen::Matrix<Scalar, Eigen::Dynamic, 1> fluence(time.rows());
  vice::functors::StarryFunctor<Scalar> func(u, r, b, tau);
  for (int i = 0; i < time.rows(); ++i) {
    Scalar lower = time(i) - 0.5*texp;
    Scalar upper = time(i) + 0.5*texp;
    fluence(i) = func.integrate(lower, upper);
  }
  return fluence;
}

template <typename Scalar>
Eigen::Matrix<Scalar, Eigen::Dynamic, 1>
transit_flux (Eigen::Ref<const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>> u, Scalar r, Scalar b, Scalar tau,
              Eigen::Ref<const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>> time)
{
  Eigen::Matrix<Scalar, Eigen::Dynamic, 1> flux(time.rows());
  vice::functors::StarryFunctor<Scalar> func(u, r, b, tau);
  for (int i = 0; i < time.rows(); ++i) {
    flux(i) = func(time(i));
  }
  return flux;
}

//
// Python interface
//
PYBIND11_MODULE(benchmark, m) {

  m.def("trapezoid_flux",            &trapezoid_flux<double>);
  m.def("transit_flux",              &transit_flux<double>);

  m.def("trapezoid_exact",           &trapezoid_exact<double>);
  m.def("transit_exact",             &transit_exact<double>);

  m.def("trapezoid_riemann",         &benchmark_trapezoid<double, unsigned, vice::integrate::riemann>);
  m.def("trapezoid_trapezoid_fixed", &benchmark_trapezoid<double, unsigned, vice::integrate::trapezoid_fixed>);
  m.def("trapezoid_simpson_fixed",   &benchmark_trapezoid<double, unsigned, vice::integrate::simpson_fixed>);
  m.def("trapezoid_trapezoid_adapt", &benchmark_trapezoid<double, double, vice::integrate::trapezoid_adapt>);
  m.def("trapezoid_simpson_adapt",   &benchmark_trapezoid<double, double, vice::integrate::simpson_adapt>);
  m.def("trapezoid_gauss",           &benchmark_trapezoid<double, double, vice::integrate::quadrature<15>>);

  m.def("transit_riemann",           &benchmark_transit<double, unsigned, vice::integrate::riemann>);
  m.def("transit_trapezoid_fixed",   &benchmark_transit<double, unsigned, vice::integrate::trapezoid_fixed>);
  m.def("transit_simpson_fixed",     &benchmark_transit<double, unsigned, vice::integrate::simpson_fixed>);
  m.def("transit_trapezoid_adapt",   &benchmark_transit<double, double, vice::integrate::trapezoid_adapt>);
  m.def("transit_simpson_adapt",     &benchmark_transit<double, double, vice::integrate::simpson_adapt>);
  m.def("transit_gauss",             &benchmark_transit<double, double, vice::integrate::quadrature<15>>);

}
