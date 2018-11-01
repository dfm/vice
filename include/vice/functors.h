#ifndef _VICE_FUNCTORS_H_
#define _VICE_FUNCTORS_H_

#include <cmath>
#include <algorithm>

#include <starry/limbdark.h>

#include "vice/integrate.h"

namespace vice {
  namespace functors {

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
    // A functor describing a starry transit
    //
    template <typename T>
    class StarryFunctor {
      private:
        long int n_eval_;

        T r_, b2_, v2_, tau_;
        limbdark::GreensLimbDark<T> L;
        utils::Vector<T> agol_c;
        T agol_norm;

      public:
        StarryFunctor (const utils::Vector<double>& u, T r, T b, T tau) : L(u.rows()), r_(r), tau_(tau) {
          b2_ = b*b;
          v2_ = 4 * (1 - b2_) / (tau*tau);

          // Convert to Agol basis
          int lmax = u.rows();
          utils::Vector<T> u_(lmax + 1);
          u_(0) = -1.0;
          u_.segment(1, lmax) = u.template cast<T>();
          agol_c = limbdark::computeC(u_);
          agol_norm = limbdark::normC(agol_c);
        }

        void reset () { n_eval_ = 0; }
        long int get_n_eval () const { return n_eval_; }
        T get_full_duration () const { return tau_; }

        T operator() (T t) {
          n_eval_++;
          T bt = sqrt(b2_ + v2_ * t*t);
          if (bt >= 1 + r_) return 1.0;
          L.compute(bt, r_);
          return L.S.dot(agol_c) * agol_norm;
        }

        T integrate (T lower, T upper) {
          auto func = [this](T t){ return this->operator()(t); };
          return vice::quadrature::integrate_adapt<15>(func, lower, upper, 1e-12, 10);
        }
    };

  }  // namespace functors
}    // namespace vice

#endif  // _VICE_FUNCTORS_H_
