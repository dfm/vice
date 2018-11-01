#ifndef _VICE_BENCHMARK_H_
#define _VICE_BENCHMARK_H_

#include <tuple>
#include <Eigen/Core>

// Timer for the benchmark.
#if defined(_MSC_VER)

//no sys/time.h in visual c++
//http://jakascorner.com/blog/2016/04/time-measurement.html
#include <chrono>
double get_timestamp() {
  using micro_s = std::chrono::microseconds;
  auto tnow = std::chrono::steady_clock::now();
  auto d_micro = std::chrono::duration_cast<micro_s>(tnow.time_since_epoch()).count();
  return double(d_micro) * 1.0e-6;
}

#else

//no std::chrono in g++ 4.8
#include <sys/time.h>
double get_timestamp() {
  struct timeval now;
  gettimeofday (&now, NULL);
  return double(now.tv_usec) * 1.0e-6 + double(now.tv_sec);
}

#endif

namespace vice {
  namespace benchmark {

    //
    // A method for benchmarking an integration method
    //
    // Returns a tuple with the computed integral, the number of calls to the
    // functor, and the total run time in seconds.
    //
    template <typename Type, typename Functor, typename Integrator, typename ParamType>
    std::tuple<Eigen::Matrix<typename Type::Scalar, Eigen::Dynamic, 1>,
               Eigen::Matrix<long int, Eigen::Dynamic, 1>, double> integrate (Functor& func, Integrator& integrator, const Eigen::DenseBase<Type>& time, const typename Type::Scalar texp, ParamType param) {
      typedef typename Type::Scalar Scalar;
      Eigen::Matrix<Scalar, Eigen::Dynamic, 1>   integral(time.rows());
      Eigen::Matrix<long int, Eigen::Dynamic, 1> n_eval(time.rows());

      double total_time = get_timestamp();
      for (int i = 0; i < time.rows(); ++i) {
        Scalar t = time(i);
        Scalar lower = t - 0.5*texp;
        Scalar upper = t + 0.5*texp;
        func.reset();
        integral(i) = integrator(func, lower, upper, param) / texp;
        n_eval(i) = func.get_n_eval();
      }
      total_time = get_timestamp() - total_time;

      return std::make_tuple(integral, n_eval, total_time);
    }

  }
}

#endif  // _VICE_BENCHMARK_H_
