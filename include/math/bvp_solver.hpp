/**
 * @file bvp_solver.hpp
 * @brief Solver for the stationary heat equation BVP.
 */

#pragma once

#include "core/config.hpp"
#include "core/types.hpp"

#include <cstddef>
#include <memory>

namespace heat {

/**
 * @brief Solves the boundary value problem using the balance method.
 *
 * This class parses the string expressions from the config using ExprTk,
 * builds the tridiagonal scheme, and solves it via the Thomas algorithm.
 */
class bvp_solver {
public:
    /**
     * @brief Result of a BVP solve on a specific grid.
     */
    struct result {
        vector x_;  ///< Grid nodes.
        vector v_;  ///< Numerical solution at the nodes.
    };

    /**
     * @brief Constructs the solver and compiles all expressions.
     * @param cfg The configuration to use.
     * @throws invalid_problem_error if any expression fails to compile.
     */
    bvp_solver(config cfg);
    ~bvp_solver();

    bvp_solver(bvp_solver&&) noexcept;
    bvp_solver& operator=(bvp_solver&&) noexcept;

    bvp_solver(const bvp_solver&) = delete;
    bvp_solver& operator=(const bvp_solver&) = delete;

    /**
     * @brief Solves the BVP on a uniform grid with n intervals.
     * @param n Number of intervals.
     * @return  The grid and the numerical solution.
     */
    [[nodiscard]] result solve(std::size_t n) const;

private:
    struct impl;
    std::unique_ptr<impl> impl_;
};

}  // namespace heat
