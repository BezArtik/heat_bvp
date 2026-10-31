/**
 * @file types.hpp
 * @brief Fundamental types for the BVP solver.
 */

#pragma once

#include <span>
#include <vector>

namespace heat {

/// Floating-point type used for all computations.
using scalar = double;

/// Vector of scalars.
using vector = std::vector<scalar>;

/// Span of scalars
using span = std::span<const scalar>;

}  // namespace heat
