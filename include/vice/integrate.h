#ifndef _VICE_INTEGRATE_H_
#define _VICE_INTEGRATE_H_

#include <cmath>
#include <tuple>
#include <algorithm>

#include <boost/math/quadrature/gauss_kronrod.hpp>

namespace vice {
  namespace integrate {

  struct riemann {

    template <typename Scalar, typename Functor>
    inline Scalar operator() (Functor& func, Scalar lower, Scalar upper, unsigned points) {
      points = std::max<unsigned>(1, points);

      Scalar dx = (upper - lower) / points;
      Scalar x = lower + 0.5 * dx;
      Scalar f = Scalar(0);
      for (unsigned i = 0; i < points; ++i) {
        f += func(x);
        x += dx;
      }

      return dx * f;
    }

  };

  struct trapezoid_adapt {

    template <typename Scalar, typename Functor>
    inline Scalar inner (Functor& func, Scalar x1, Scalar x2, Scalar y1, Scalar y2, Scalar tol, unsigned max_depth, unsigned depth) {
      Scalar x0 = 0.5 * (x1 + x2);
      Scalar val = func(x0);
      Scalar dx = (x2 - x1);
      Scalar pred = 0.5 * dx * (y1 + y2);
      Scalar integral = 0.25 * dx * (y1 + 2*val + y2);

      if (depth < max_depth && std::abs(pred - integral) > tol) {
        integral = inner<Scalar, Functor>(func, x1, x0, y1, val, tol, max_depth, depth+1);
        integral += inner<Scalar, Functor>(func, x0, x2, val, y2, tol, max_depth, depth+1);
      }

      return integral;
    }

    template <typename Scalar, typename Functor>
    inline Scalar operator() (Functor& func, Scalar lower, Scalar upper, Scalar tol, unsigned max_depth=15) {
      return inner<Scalar, Functor> (func, lower, upper, func(lower), func(upper), tol, max_depth, 0);
    }

  };

  struct trapezoid_fixed {

    template <typename Scalar, typename Functor>
    inline Scalar operator() (Functor& func, Scalar lower, Scalar upper, unsigned points) {
      points = std::max<unsigned>(2, points);

      Scalar dx = (upper - lower) / (points - 1);
      Scalar x = lower + dx;
      Scalar f = 0.5 * func(lower);
      for (unsigned i = 1; i < points - 1; ++i) {
        f += func(x);
        x += dx;
      }
      f += 0.5 * func(upper);

      return dx * f;
    }

  };

  struct simpson_adapt {

    template <typename Scalar, typename Functor>
    inline Scalar inner (Functor& func, Scalar x0, Scalar dx, Scalar ym, Scalar y0, Scalar yp, Scalar tol, unsigned max_depth, unsigned depth) {
      Scalar x_m = x0 - 0.5*dx;
      Scalar val_m = func(x_m);
      Scalar int_m = dx * (4*val_m + ym + y0) / 6;

      Scalar x_p = x0 + 0.5*dx;
      Scalar val_p = func(x_p);
      Scalar int_p = dx * (4*val_p + yp + y0) / 6;

      Scalar pred = dx * (4*y0 + ym + yp) / 3;

      if (depth < max_depth && std::abs(pred - (int_m + int_p)) > 15 * tol) {
        int_m = inner<Scalar, Functor>(func, x_m, 0.5*dx, ym, val_m, y0, tol, max_depth, depth+1);
        int_p = inner<Scalar, Functor>(func, x_p, 0.5*dx, y0, val_p, yp, tol, max_depth, depth+1);
      }

      return int_m + int_p;
    }

    template <typename Scalar, typename Functor>
    inline Scalar operator() (Functor& func, Scalar lower, Scalar upper, Scalar tol, unsigned max_depth=15) {
      Scalar x0 = 0.5 * (upper + lower);
      Scalar dx = 0.5*(upper - lower);
      Scalar ym = func(lower);
      Scalar y0 = func(x0);
      Scalar yp = func(upper);

      Scalar int_trap = 0.5 * dx * (ym + 2 * y0 + yp);
      Scalar int_simp = dx * (ym + 4 * y0 + yp) / 3;

      if (std::abs(int_trap - int_simp) < tol) {
        return int_simp;
      }

      return inner<Scalar, Functor> (func, x0, dx, ym, y0, yp, tol, max_depth, 0);
    }
  };

  struct simpson_fixed {

    template <typename Scalar, typename Functor>
    inline Scalar operator() (Functor& func, Scalar lower, Scalar upper, unsigned points) {
      points = std::max<unsigned>(3, points + (points + 1) % 2);  // points must be odd

      Scalar dx = (upper - lower) / (points - 1);
      Scalar x = lower + dx;
      Scalar f = func(lower);
      for (unsigned i = 1; i < points - 1; ++i) {
        f += (2 + 2 * (i % 2)) * func(x);
        x += dx;
      }
      f += func(upper);

      return dx * f / 3;
    }

  };

  template <unsigned Points>
  struct quadrature {

    template <typename Scalar, typename Functor>
    inline Scalar operator() (Functor& func, Scalar lower, Scalar upper, Scalar tol, unsigned max_depth=15) {
      auto f = [&func](Scalar t) { return func(t); };
      return boost::math::quadrature::gauss_kronrod<Scalar, Points>::integrate(f, lower, upper, max_depth, tol);
    }

  };

  }  // namespace integrate
}    // namespace vice

#endif  // _VICE_INTEGRATE_H_
