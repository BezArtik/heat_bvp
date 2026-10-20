/**
 * @file config.hpp
 * @brief Configuration structure for the BVP solver.
 */

#pragma once

#include "core/types.hpp"

#include <cstddef>
#include <string>

namespace heat {

/**
 * @brief Expressions for a specific side of the discontinuity.
 */
struct side_functions {
    std::string k_;  ///< Thermal conductivity k(x).
    std::string q_;  ///< Heat source q(x).
    std::string f_;  ///< Right-hand side f(x).
};

/**
 * @brief Full configuration for a single BVP run.
 */
struct config {
    // Problem definition
    scalar a_ = 0.0;
    scalar b_ = 1.0;
    scalar xi_ = 0.5;
    scalar mu1_ = 0.0;
    scalar mu2_ = 0.0;
    std::string boundary_type_ = "dirichlet";
    std::string analytical_;

    side_functions left_;
    side_functions right_;

    // Solver settings
    std::size_t n_ = 10;
    scalar eps_ = 5e-7;
    std::size_t max_n_ = 100000;

    // Output
    std::string output_path_;
    std::size_t print_every_ = 1;
};

}  // namespace heat
