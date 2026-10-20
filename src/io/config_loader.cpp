/**
 * @file config_loader.cpp
 * @brief Implementation of the TOML configuration loader.
 */

#include "io/config_loader.hpp"

#include "core/errors.hpp"

#include <format>
#include <string>
#include <string_view>
#include <toml++/toml.hpp>

namespace heat {

namespace {

[[noreturn]] void missing_field(std::string_view path) {
    throw invalid_problem_error{std::format("config: missing required field '{}'", path)};
}

[[noreturn]] void wrong_type(std::string_view path, std::string_view expected) {
    throw invalid_problem_error{std::format("config: field '{}' must be {}", path, expected)};
}

// --- Typed accessors ---

scalar get_float(const toml::table& t, std::string_view path) {
    auto view = t.at_path(path);
    auto node = view.as_floating_point();
    if (!node) {
        if (view) wrong_type(path, "a floating-point number");
        missing_field(path);
    }
    return node->get();
}

std::size_t get_size(const toml::table& t, std::string_view path) {
    auto view = t.at_path(path);
    auto node = view.as_integer();
    if (!node) {
        if (view) wrong_type(path, "an integer");
        missing_field(path);
    }
    const auto val = node->get();
    if (val < 0) throw invalid_problem_error{std::format("config: '{}' must be non-negative", path)};

    return val;
}

std::string get_string(const toml::table& t, std::string_view path) {
    auto view = t.at_path(path);
    auto node = view.as_string();
    if (!node) {
        if (view) wrong_type(path, "a string");
        missing_field(path);
    }
    return node->get();
}

}  // namespace

config load_config(std::string_view path) {
    auto tbl = toml::table{};
    try {
        tbl = toml::parse_file(path);
    } catch (const toml::parse_error& e) {
        throw invalid_problem_error{std::format("config: failed to parse '{}': {}", path, e.description())};
    }

    auto cfg = config{};

    // --- Problem section ---
    cfg.a_ = get_float(tbl, "problem.a");
    cfg.b_ = get_float(tbl, "problem.b");
    cfg.xi_ = get_float(tbl, "problem.xi");
    cfg.mu1_ = get_float(tbl, "problem.mu1");
    cfg.mu2_ = get_float(tbl, "problem.mu2");

    cfg.boundary_type_ = get_string(tbl, "problem.boundary_type");
    cfg.analytical_ = get_string(tbl, "problem.analytical");

    cfg.left_.k_ = get_string(tbl, "problem.left.k");
    cfg.left_.q_ = get_string(tbl, "problem.left.q");
    cfg.left_.f_ = get_string(tbl, "problem.left.f");

    cfg.right_.k_ = get_string(tbl, "problem.right.k");
    cfg.right_.q_ = get_string(tbl, "problem.right.q");
    cfg.right_.f_ = get_string(tbl, "problem.right.f");

    // --- Solver section ---
    cfg.n_ = get_size(tbl, "solver.n");
    cfg.eps_ = get_float(tbl, "solver.eps");
    cfg.max_n_ = get_size(tbl, "solver.max_n");

    // --- Output section ---
    cfg.output_path_ = get_string(tbl, "output.path");
    cfg.print_every_ = get_size(tbl, "output.print_every");

    if (cfg.print_every_ == 0) throw invalid_problem_error{"config: 'output.print_every' must be at least 1"};

    return cfg;
}

}  // namespace heat
