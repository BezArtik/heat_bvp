#include "app/app.hpp"

#include "core/errors.hpp"
#include "exprtk.hpp"
#include "math/bvp_solver.hpp"

#include <algorithm>
#include <cmath>
#include <format>
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
    scalar operator()(scalar x) const {
        x_ = x;
        return expr_.value();
    }

private:
    mutable scalar x_{};
    exprtk::symbol_table<scalar> symtab_;
    exprtk::expression<scalar> expr_;
};

solve_result run_test(const config& cfg, bvp_solver& solver) {
    auto u_exact = expr_scalar{cfg.analytical_, cfg.xi_};
    auto iter = std::size_t{}, n = cfg.n_;

    while (true) {
        ++iter;
        if (n > cfg.max_n_)
            throw solver_error{std::format("test: max_n ({}) reached, eps = {} not achieved", cfg.max_n_, cfg.eps_)};

        auto sol = solver.solve(n);

        using pair = std::pair<scalar, scalar>;
        auto err = std::views::zip(sol.x_, sol.v_) | std::views::transform([&](auto&& p) -> pair {
                       auto&& [x, v] = p;
                       return {std::abs(v - u_exact(x)), x};
                   });
        auto&& [eps, eps_at] = *std::ranges::max_element(err, {}, &pair::first);

        if (eps <= cfg.eps_ || cfg.eps_ < 0.0) {
            auto reference = vector(sol.v_.size());

            std::ranges::transform(sol.x_, reference.begin(), [&](auto&& x) { return u_exact(x); });
            auto diff = std::views::zip(sol.v_, reference) | std::views::transform([](auto&& p) {
                            auto&& [v, ref] = p;
                            return v - ref;
                        });
            auto error = vector(diff.begin(), diff.end());

            return {.x_ = std::move(sol.x_),
                    .v_ = std::move(sol.v_),
                    .reference_ = std::move(reference),
                    .error_ = std::move(error),
                    .final_n_ = n,
                    .iterations_ = iter,
                    .eps_achieved_ = eps,
                    .eps_at_ = eps_at,
                    .is_test_ = true};
        }
        n *= 2;
    }
}

solve_result run_main(const config& cfg, bvp_solver& solver) {
    auto iter = std::size_t{}, n = cfg.n_;
    auto v_coarse = solver.solve(n);
    auto v_fine = solver.solve(2 * n);

    while (true) {
        ++iter;

        auto v_fine_strided = v_fine.v_ | std::views::stride(2);

        using pair = std::pair<scalar, scalar>;
        auto err =
            std::views::zip(v_coarse.x_, v_coarse.v_, v_fine_strided) | std::views::transform([](auto&& p) -> pair {
                auto&& [x, v, v2] = p;
                return {std::abs(v - v2), x};
            });

        auto&& [eps, eps_at] = *std::ranges::max_element(err, {}, &pair::first);

        if (eps <= cfg.eps_ || cfg.eps_ < 0.0) {
            auto diff = std::views::zip(v_coarse.v_, v_fine_strided) | std::views::transform([](auto&& p) {
                            auto&& [v, ref] = p;
                            return v - ref;
                        });
            auto reference = vector{v_fine_strided.begin(), v_fine_strided.end()};
            auto error = vector{diff.begin(), diff.end()};

            return {.x_ = std::move(v_coarse.x_),
                    .v_ = std::move(v_coarse.v_),
                    .reference_ = std::move(reference),
                    .error_ = std::move(error),
                    .final_n_ = n,
                    .iterations_ = iter,
                    .eps_achieved_ = eps,
                    .eps_at_ = eps_at,
                    .is_test_ = false};
        }

        n *= 2;
        if (2 * n > cfg.max_n_) throw solver_error{std::format("main: max_n reached, eps = {} not achieved", cfg.eps_)};

        v_coarse = std::move(v_fine);
        v_fine = solver.solve(2 * n);
    }
}

}  // namespace

solve_result solve_with_refinement(const config& cfg) {
    auto solver = bvp_solver{cfg};
    return !cfg.analytical_.empty() ? run_test(cfg, solver) : run_main(cfg, solver);
}

}  // namespace heat
