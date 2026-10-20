/**
 * @file config_loader.hpp
 * @brief Loads configuration from TOML files.
 */

#pragma once

#include "core/config.hpp"

#include <string_view>

namespace heat {

/**
 * @brief Loads and validates a configuration from a TOML file.
 *
 * @param path Path to the TOML file.
 * @return Parsed configuration.
 *
 * @throws invalid_problem_error if the file cannot be read,
 *         contains syntax errors, or misses required fields.
 */
[[nodiscard]] config load_config(std::string_view path);

}  // namespace heat
