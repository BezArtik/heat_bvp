/**
 * @file scheme.hpp
 * @brief Construction of the finite volume (balance method) scheme.
 */

#pragma once

#include "core/types.hpp"
#include "math/grid.hpp"

#include <concepts>
#include <cstddef>
#include <functional>
#include <ranges>

namespace heat {

/**
 * @brief Represents a tridiagonal linear system.
 *
 * Equations are of the form:
 * lower[i] * v[i-1] + main[i] * v[i] + upper[i] * v[i+1] = rhs[i]
 */
struct tridiagonal_system {
    vector lower_;  ///< Lower diagonal (size n). lower[0] is unused for Dirichlet.
    vector main_;   ///< Main diagonal (size n+1).
    vector upper_;  ///< Upper diagonal (size n). upper[n-1] is unused for Dirichlet.
    vector rhs_;    ///< Right-hand side (size n+1).
};

template <typename F>
concept integrate_callable = std::invocable<F, scalar>;

/**
 * @brief Integrates a function using the composite Simpson's rule.
 *
 * @tparam F Callable type: scalar(scalar).
 * @param f      The function to integrate.
 * @param a      Left bound.
 * @param b      Right bound.
 * @param n_sub  Number of sub-intervals (must be even). Default is 2.
 * @return       Integral of f from a to b.
 */
template <integrate_callable F>
    requires std::convertible_to<std::invoke_result_t<F, scalar>, scalar>
scalar integrate_simpson(F f, scalar a, scalar b, std::size_t n_sub = 2) {
    if (a == b) return scalar{};

    if (n_sub % 2 != 0) ++n_sub;

    const auto h = (b - a) / static_cast<scalar>(n_sub);

    auto sum = std::invoke(f, a) + std::invoke(f, b);

    for (const auto i : std::views::iota(1ul, n_sub)) {
        const auto val = std::invoke(f, a + static_cast<scalar>(i) * h);
        sum += (i % 2 == 0) ? 2.0 * val : 4.0 * val;
    }

    return sum * h / 3.0;
}

/**
 * @brief Integrates a piecewise function on [a, b].
 *
 * If the interval [a, b] contains the breakpoint xi, the integral is
 * split into two parts: [a, xi] using f_left, and [xi, b] using f_right.
 *
 * @tparam FLeft  Callable for the left side.
 * @tparam FRight Callable for the right side.
 */
template <integrate_callable FLeft, integrate_callable FRight>
scalar integrate(const grid& g, FLeft f_left, FRight f_right, scalar a, scalar b) {
    const scalar xi = g.xi();
    if (b <= xi) {
        return integrate_simpson(f_left, a, b);
    } else if (a >= xi) {
        return integrate_simpson(f_right, a, b);
    } else {
        return integrate_simpson(f_left, a, xi) + integrate_simpson(f_right, xi, b);
    }
}

/**
 * @brief Assembles the tridiagonal system for the BVP.
 *
 * The problem is:
 *   (k(x)u')' - q(x)u = -f(x)
 *
 * with Dirichlet boundary conditions u(a) = mu1, u(b) = mu2.
 *
 * @tparam K1,K2 Callables for k(x) on the left and right of xi.
 * @tparam Q1,Q2 Callables for q(x) on the left and right of xi.
 * @tparam F1,F2 Callables for f(x) on the left and right of xi.
 *
 * @param g   The uniform grid.
 * @param k1  k(x) for x < xi.
 * @param k2  k(x) for x > xi.
 * @param q1  q(x) for x < xi.
 * @param q2  q(x) for x > xi.
 * @param f1  f(x) for x < xi.
 * @param f2  f(x) for x > xi.
 * @param mu1 Value at x = a.
 * @param mu2 Value at x = b.
 * @return    The assembled tridiagonal system.
 */
template <integrate_callable K1, integrate_callable K2, integrate_callable Q1, integrate_callable Q2,
          integrate_callable F1, integrate_callable F2>
auto assemble_scheme(const grid& g, K1 k1, K2 k2, Q1 q1, Q2 q2, F1 f1, F2 f2, scalar mu1, scalar mu2) {
    const auto n = g.n();
    const auto h = g.h();
    const auto h2 = h * h;

    auto sys =
        tridiagonal_system{.lower_ = vector(n), .main_ = vector(n + 1), .upper_ = vector(n), .rhs_ = vector(n + 1)};

    // --- Boundary conditions ---
    // Left node (i = 0)
    sys.main_[0] = 1.0;
    sys.rhs_[0] = mu1;

    // Right node (i = n)
    sys.main_[n] = 1.0;
    sys.rhs_[n] = mu2;

    // --- Precompute the inverse of k(x) ---
    auto inv_k1 = [&k1](scalar x) { return 1.0 / std::invoke(k1, x); };
    auto inv_k2 = [&k2](scalar x) { return 1.0 / std::invoke(k2, x); };

    // --- Compute a_i coefficients for i = 1..n ---
    // a_i = h / ( \int_{x_{i-1}}^{x_i} dx / k(x) )
    auto a = vector(n + 1);
    for (const auto i : std::views::iota(1ul, n + 1)) a[i] = h / integrate(g, inv_k1, inv_k2, g.x(i - 1), g.x(i));

    // --- Assemble internal equations (i = 1..n-1) ---
    for (const auto i : std::views::iota(1ul, n)) {
        // d_i = (1/h) * \int_{x_{i-0.5}}^{x_{i+0.5}} q(x) dx
        const auto d_i = integrate(g, q1, q2, g.x_half(i - 1), g.x_half(i)) / h;

        // phi_i = (1/h) * \int_{x_{i-0.5}}^{x_{i+0.5}} f(x) dx
        const auto phi_i = integrate(g, f1, f2, g.x_half(i - 1), g.x_half(i)) / h;

        // The equation:
        // (a_i/h^2) v_{i-1} - ((a_i+a_{i+1})/h^2 + d_i) v_i + (a_{i+1}/h^2) v_{i+1} = -phi_i
        sys.lower_[i - 1] = a[i] / h2;
        sys.main_[i] = -(a[i] + a[i + 1]) / h2 - d_i;
        sys.upper_[i] = a[i + 1] / h2;
        sys.rhs_[i] = -phi_i;
    }

    return sys;
}

}  // namespace heat
