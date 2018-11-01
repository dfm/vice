#ifndef _VICE_BENCHMARK_H_
#define _VICE_BENCHMARK_H_

#include <cmath>
#include <algorithm>

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

//
// A functor describing a trapezoidal transit
//
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


//
// A method for benchmarking an integration method
//
// Returns a tuple with the computed integral, the number of calls to the
// functor, and the total run time in seconds.
//
template <typename Scalar, typename Functor, typename Integrator, typename ParamType>
std::tuple<std::vector<Scalar>, long int, double> benchmark_integrate (Functor& func, Integrator& integrator, const std::vector<Scalar>& time, const Scalar texp, ParamType param) {
  func.reset();
  std::vector<Scalar> fluence(time.size());

  double total_time = get_timestamp();
  for (size_t i = 0; i < time.size(); ++i) {
    Scalar t = time[i];
    Scalar lower = t - 0.5*texp;
    Scalar upper = t + 0.5*texp;
    fluence[i] = integrator(func, lower, upper, param) / texp;
  }
  total_time = get_timestamp() - total_time;

  return std::make_tuple(fluence, func.get_n_eval(), total_time);
}

#endif  // _VICE_BENCHMARK_H_
