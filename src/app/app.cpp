#include "app/app.hpp"

#include "core/errors.hpp"
#include "exprtk.hpp"
#include "io/config_loader.hpp"
#include "math/bvp_solver.hpp"

#include <algorithm>
#include <cmath>
#include <format>
#include <print>
#include <ranges>
#include <utility>

namespace heat {

namespace {

// Helper to evaluate an analytical expression at scalar points.
class expr_scalar {
public:
    expr_scalar(const std::string& src, scalar xi) {
        symtab_.add_variable("x", x_);
        symtab_.add_constant("xi", xi);
        expr_.register_symbol_table(symtab_);
        auto parser = exprtk::parser<scalar>{};
        if (!parser.compile(src, expr_))
            throw invalid_problem_error{std::format("failed to compile analytical expression: {}", parser.error())};
    }
    auto operator()(scalar x) const {
        x_ = x;
        return expr_();
    }

private:
    mutable scalar x_{};
    exprtk::symbol_table<scalar> symtab_;
    exprtk::expression<scalar> expr_;
};

// Result of solving a BVP with automatic grid refinement.
struct solve_result {
    vector x_;
    vector v_;
    vector reference_;
    vector error_;

    std::size_t final_n_;
    std::size_t iterations_;
    scalar eps_achieved_;
    scalar eps_at_;
    bool is_test_;
};

// (max|v - ref|, x_at_max) for parallel ranges x, v and ref.
using eps_pair = std::pair<scalar, scalar>;

namespace rng = std::ranges;
namespace view = rng::views;

auto find_max_error(const bvp_solver::result& res, span ref) {
    auto err = view::zip(res.x_, res.v_, ref) | view::transform([](auto&& p) -> eps_pair {
                   auto&& [x, v, r] = p;
                   return {std::abs(v - r), x};
               });
    return *rng::max_element(err, {}, &eps_pair::first);
}

// Packs a fully-solved state into the result struct.
auto finalize(bvp_solver::result res, vector ref, std::size_t n, std::size_t iter, eps_pair eps_info, bool is_test) {
    auto error = view::zip(res.v_, ref) | view::transform([](auto&& p) {
                     auto&& [a, b] = p;
                     return a - b;
                 }) |
                 rng::to<vector>();

    return solve_result{.x_ = std::move(res.x_),
                        .v_ = std::move(res.v_),
                        .reference_ = std::move(ref),
                        .error_ = std::move(error),
                        .final_n_ = n,
                        .iterations_ = iter,
                        .eps_achieved_ = eps_info.first,
                        .eps_at_ = eps_info.second,
                        .is_test_ = is_test};
}

auto run_test(const config& cfg, const bvp_solver& solver) {
    const auto u_exact = expr_scalar{cfg.analytical_, cfg.xi_};

    for (auto iter = std::size_t{1}, n = cfg.n_;; ++iter, n *= 2) {
        if (n > cfg.max_n_)
            throw solver_error{std::format("test: max_n ({}) reached, eps = {} not achieved", cfg.max_n_, cfg.eps_)};

        auto sol = solver.solve(n);
        auto ref = sol.x_ | view::transform([&](auto x) { return u_exact(x); }) | rng::to<vector>();
        const auto eps = find_max_error(sol, ref);

        if (eps.first <= cfg.eps_ || cfg.eps_ < 0.0)
            return finalize(std::move(sol), std::move(ref), n, iter, eps, true);
    }
}

auto run_main(const config& cfg, const bvp_solver& solver) {
    auto n = cfg.n_;
    auto v_coarse = solver.solve(n);
    auto v_fine = solver.solve(2 * n);

    for (auto iter = std::size_t{1};; ++iter) {
        auto fine_strided = v_fine.v_ | view::stride(2) | rng::to<vector>();
        const auto eps = find_max_error(v_coarse, fine_strided);

        if (eps.first <= cfg.eps_ || cfg.eps_ < 0.0)
            return finalize(std::move(v_coarse), std::move(fine_strided), n, iter, eps, false);

        n *= 2;
        if (2 * n > cfg.max_n_) throw solver_error{std::format("main: max_n reached, eps = {} not achieved", cfg.eps_)};

        v_coarse = std::move(v_fine);
        v_fine = solver.solve(2 * n);
    }
}

}  // namespace

int run(int argc, char** argv) {
    if (argc != 2) {
        std::print(stderr, "Usage: {} <config.toml>\n", argv[0]);
        return 1;
    }

    try {
        const auto cfg = load_config(argv[1]);
        const auto solver = bvp_solver{cfg};
        const auto res = !cfg.analytical_.empty() ? run_test(cfg, solver) : run_main(cfg, solver);

        auto file = std::ofstream{cfg.output_path_};
        if (!file) throw solver_error{std::format("cannot open '{}'", cfg.output_path_)};

        std::print(file, "{}\n", res.is_test_ ? "i,x,u,v,err" : "i,x,v,v2,err");
        for (const auto i : view::iota(0ul, res.x_.size()))
            std::print(file, "{},{},{},{},{}\n", i, res.x_[i], res.v_[i], res.reference_[i], res.error_[i]);

        std::print(stderr,
                   "n = {}\n"
                   "iterations = {}\n"
                   "{} = {}\n"
                   "max at x = {}\n",
                   res.final_n_, res.iterations_, (res.is_test_ ? "eps1" : "eps2"), res.eps_achieved_, res.eps_at_);

        return 0;
    } catch (const heat::solver_error& e) {
        std::print(stderr, "Error: {}\n", e.what());
        return 1;
    } catch (const std::exception& e) {
        std::print(stderr, "Unexpected: {}\n", e.what());
        return 2;
    }
}

}  // namespace heat
