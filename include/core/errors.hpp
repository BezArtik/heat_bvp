/**
 * @file errors.hpp
 * @brief Exception hierarchy for the BVP solver.
 */

#pragma once

#include <stdexcept>

namespace heat {

/**
 * @brief Base class for all solver errors.
 */
struct solver_error : std::runtime_error {
    using std::runtime_error::runtime_error;
};

/**
 * @brief Thrown when the problem definition is invalid.
 */
struct invalid_problem_error : solver_error {
    using solver_error::solver_error;
};

}  // namespace heat
