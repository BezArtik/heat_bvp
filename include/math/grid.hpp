/**
 * @file grid.hpp
 * @brief Uniform grid definition.
 */

#pragma once

#include "core/types.hpp"

#include <cmath>
#include <cstddef>
#include <limits>

namespace heat {

/**
 * @brief Represents a uniform grid on the interval [a, b].
 *
 * The grid consists of n intervals and n+1 nodes.
 * It also stores the breakpoint xi, which may or may not
 * coincide with one of the nodes.
 */
class grid {
public:
    /**
     * @brief Constructs a uniform grid.
     * @param a  Left boundary.
     * @param b  Right boundary.
     * @param n  Number of intervals.
     * @param xi Breakpoint of the coefficients.
     */
    grid(scalar a, scalar b, std::size_t n, scalar xi) noexcept
        : a_{a}, b_{b}, n_{n}, xi_{xi}, h_{(b - a) / static_cast<scalar>(n)} {
        const auto diff = (xi_ - a_) / h_;
        is_xi_at_node_ = (std::abs(diff - std::round(diff)) < std::numeric_limits<scalar>::epsilon());
        xi_index_ = is_xi_at_node_ ? static_cast<std::size_t>(std::round(diff)) : 0;
    }

    /// Returns the number of intervals.
    [[nodiscard]] std::size_t n() const noexcept { return n_; }

    /// Returns the number of nodes (n + 1).
    [[nodiscard]] std::size_t size() const noexcept { return n_ + 1; }

    /// Returns the step size h.
    [[nodiscard]] scalar h() const noexcept { return h_; }

    /// Returns the left boundary.
    [[nodiscard]] scalar a() const noexcept { return a_; }

    /// Returns the right boundary.
    [[nodiscard]] scalar b() const noexcept { return b_; }

    /// Returns the breakpoint xi.
    [[nodiscard]] scalar xi() const noexcept { return xi_; }

    /**
     * @brief Returns the coordinate of the i-th node.
     * @param i Node index (0 <= i <= n).
     */
    [[nodiscard]] scalar x(std::size_t i) const noexcept { return a_ + static_cast<scalar>(i) * h_; }

    /**
     * @brief Returns the coordinate of the i-th half-node.
     * @param i Half-node index (0 <= i < n).
     */
    [[nodiscard]] scalar x_half(std::size_t i) const noexcept { return a_ + (static_cast<scalar>(i) + 0.5) * h_; }

    /**
     * @brief Checks if the breakpoint xi exactly matches a node.
     *
     * Due to floating-point issues, we check the distance modulo h.
     */
    [[nodiscard]] bool is_xi_at_node() const noexcept { return is_xi_at_node_; }

    /**
     * @brief Returns the node index where xi is located, if is_xi_at_node() is
     * true.
     */
    [[nodiscard]] std::size_t xi_index() const noexcept { return xi_index_; }

private:
    scalar a_;
    scalar b_;
    std::size_t n_;
    scalar xi_;
    scalar h_;
    bool is_xi_at_node_;
    std::size_t xi_index_;
};

}  // namespace heat
