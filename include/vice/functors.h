#ifndef _VICE_FUNCTORS_H_
#define _VICE_FUNCTORS_H_

#include <cmath>
#include <vector>
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

        std::vector<T> cuts_;

      public:

        TrapFunctor (T depth, T duration, T ingress)
          : depth_(depth), duration_(duration), ingress_(ingress), n_eval_(0)
        {
          cuts_.resize(4);
          cuts_[0] = -0.5*duration_ - ingress_;
          cuts_[1] = -0.5*duration_;
          cuts_[2] =  0.5*duration_;
          cuts_[3] =  0.5*duration_ + ingress_;
        }

        void reset () { n_eval_ = 0; }
        long int get_n_eval () const { return n_eval_; }
        T get_full_duration () const { return duration_ + ingress_; }

        int get_cuts (T lower, T upper, std::vector<T>& cuts) {
          int j = 1;
          cuts.resize(6);
          cuts[0] = lower;
          for (int i = 0; i < 4; ++i) {
            if (lower < cuts_[i] && cuts_[i] < upper) {
              cuts[j] = cuts_[i];
              ++j;
            }
          }
          cuts[j] = upper;
          cuts.resize(j + 1);
          return j + 1;
        }

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
        std::vector<T> cuts_;

        T opr2_, r_, b2_, v2_, tau_;
        limbdark::GreensLimbDark<T> L;
        utils::Vector<T> agol_c;
        T agol_norm;

      public:
        StarryFunctor (const utils::Vector<T>& u, T r, T b, T tau) : r_(r), tau_(tau), L(u.rows())  {
          b2_ = b*b;
          opr2_ = (1 + r_) * (1 + r_);
          v2_ = 4 * (opr2_ - b2_) / (tau*tau);

          cuts_.resize(4);
          cuts_[0] = -sqrt((opr2_ - b2_) / v2_);
          cuts_[1] = -sqrt(((1 - r_) * (1 - r_) - b2_) / v2_);
          cuts_[2] =  sqrt(((1 - r_) * (1 - r_) - b2_) / v2_);
          cuts_[3] =  sqrt((opr2_ - b2_) / v2_);
          std::sort(cuts_.begin(), cuts_.end());

          // Convert to Agol basis
          int lmax = u.rows();
          utils::Vector<T> u_(lmax + 1);
          u_(0) = -1.0;
          u_.segment(1, lmax) = u;
          agol_c = limbdark::computeC(u_);
          agol_norm = limbdark::normC(agol_c);
        }

        void reset () { n_eval_ = 0; }
        long int get_n_eval () const { return n_eval_; }
        T get_full_duration () const { return tau_; }

        int get_cuts (T lower, T upper, std::vector<T>& cuts) {
          int j = 1;
          cuts.resize(6);
          cuts[0] = lower;
          for (int i = 0; i < 4; ++i) {
            if (lower < cuts_[i] && cuts_[i] < upper) {
              cuts[j] = cuts_[i];
              ++j;
            }
          }
          cuts[j] = upper;
          cuts.resize(j + 1);
          return j + 1;
        }

        T operator() (T t) {
          n_eval_++;
          T bt2 = b2_ + v2_ * t*t;
          if (bt2 < 0.0 || bt2 >= opr2_) return 1.0;
          L.compute(sqrt(bt2), r_);
          return L.S.dot(agol_c) * agol_norm;
        }

        T integrate (T lower, T upper, T tol=1e-10, unsigned max_depth=15) {
          vice::integrate::quadrature<15> integrator;
          return integrator(*this, lower, upper, tol, max_depth);
        }
    };

  }  // namespace functors
}    // namespace vice

#endif  // _VICE_FUNCTORS_H_
