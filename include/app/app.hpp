/**
 * @file app.hpp
 * @brief Orchestration for the BVP solver.
 */

#pragma once

#include "core/config.hpp"
#include "core/types.hpp"

namespace heat {

/**
 * @brief Result of solving a BVP with automatic grid refinement.
 */
struct solve_result {
    vector x_;          ///< Grid nodes.
    vector v_;          ///< Numerical solution on the final grid.
    vector reference_;  ///< Reference values (u for test, v2 for main).
    vector error_;      ///< Signed difference (v - reference).

    std::size_t final_n_;     ///< Final grid intervals.
    std::size_t iterations_;  ///< Number of refinements.
    scalar eps_achieved_;     ///< Achieved accuracy (eps1 or eps2).
    scalar eps_at_;           ///< x coordinate of maximum error.
    bool is_test_;            ///< true if test mode (vs analytical).
};

/**
 * @brief Solves the BVP with automatic grid refinement.
 *
 * If cfg.analytical is not empty, runs in test mode (eps1 criterion).
 * Otherwise, runs in main mode (eps2 criterion with half-step).
 *
 * @param cfg Problem configuration.
 * @return Result containing the final grid, solution, and statistics.
 * @throws solver_error if max_n is reached without achieving eps.
 */
[[nodiscard]] solve_result solve_with_refinement(const config& cfg);

}  // namespace heat
