#include "math/bvp_solver.hpp"

#include "core/errors.hpp"
#include "exprtk.hpp"
#include "math/grid.hpp"
#include "math/scheme.hpp"

#include <format>
#include <utility>

namespace heat {

namespace {

auto solve_thomas(const tridiagonal_system& sys) {
    const auto n = sys.main_.size();
    if (n == 0) return vector{};
    if (sys.lower_.size() != n - 1 || sys.upper_.size() != n - 1 || sys.rhs_.size() != n)
        throw std::invalid_argument("thomas: inconsistent matrix sizes");

    auto c_prime = vector(n);
    auto y = vector(n);

    c_prime[0] = sys.upper_[0] / sys.main_[0];
    y[0] = sys.rhs_[0] / sys.main_[0];

    for (const auto i : std::views::iota(1ul, n - 1)) {
        const auto denom = sys.main_[i] - sys.lower_[i - 1] * c_prime[i - 1];
        c_prime[i] = sys.upper_[i] / denom;
        y[i] = (sys.rhs_[i] - sys.lower_[i - 1] * y[i - 1]) / denom;
    }

    y[n - 1] =
        (sys.rhs_[n - 1] - sys.lower_[n - 2] * y[n - 2]) / (sys.main_[n - 1] - sys.lower_[n - 2] * c_prime[n - 2]);

    for (std::size_t i = n - 1; i > 0; --i) y[i - 1] -= c_prime[i - 1] * y[i];

    return y;
}

}  // namespace

struct bvp_solver::impl {
    config cfg_;

    mutable scalar x_var_{};
    exprtk::symbol_table<scalar> symtab_;

    exprtk::expression<scalar> k1_expr_, k2_expr_, q1_expr_, q2_expr_, f1_expr_, f2_expr_;

    impl(config c) : cfg_{std::move(c)} {
        symtab_.add_variable("x", x_var_);
        symtab_.add_constant("xi", cfg_.xi_);

        auto compile = [this](auto&& str, auto&& expr, auto&& name) {
            expr.register_symbol_table(symtab_);
            auto parser = exprtk::parser<scalar>{};
            if (!parser.compile(str, expr))
                throw invalid_problem_error{std::format("failed to compile {}: {}", name, parser.error())};
        };

        compile(cfg_.left_.k_, k1_expr_, "left.k");
        compile(cfg_.left_.q_, q1_expr_, "left.q");
        compile(cfg_.left_.f_, f1_expr_, "left.f");

        compile(cfg_.right_.k_, k2_expr_, "right.k");
        compile(cfg_.right_.q_, q2_expr_, "right.q");
        compile(cfg_.right_.f_, f2_expr_, "right.f");
    }

    auto eval(exprtk::expression<scalar>& expr, scalar x) const {
        x_var_ = x;
        return expr();
    }
};

bvp_solver::bvp_solver(config cfg) : impl_{std::make_unique<impl>(std::move(cfg))} {}

bvp_solver::~bvp_solver() = default;
bvp_solver::bvp_solver(bvp_solver&&) noexcept = default;
bvp_solver& bvp_solver::operator=(bvp_solver&&) noexcept = default;

bvp_solver::result bvp_solver::solve(std::size_t n) const {
    auto&& impl = *impl_;
    auto&& cfg = impl.cfg_;

    if (cfg.boundary_type_ != "dirichlet")
        throw invalid_problem_error{std::format("unsupported boundary type: '{}'", cfg.boundary_type_)};

    const auto g = grid{cfg.a_, cfg.b_, n, cfg.xi_};

    auto k1 = [&impl](scalar x) { return impl.eval(impl.k1_expr_, x); };
    auto k2 = [&impl](scalar x) { return impl.eval(impl.k2_expr_, x); };
    auto q1 = [&impl](scalar x) { return impl.eval(impl.q1_expr_, x); };
    auto q2 = [&impl](scalar x) { return impl.eval(impl.q2_expr_, x); };
    auto f1 = [&impl](scalar x) { return impl.eval(impl.f1_expr_, x); };
    auto f2 = [&impl](scalar x) { return impl.eval(impl.f2_expr_, x); };

    const auto sys = assemble_scheme(g, k1, k2, q1, q2, f1, f2, cfg.mu1_, cfg.mu2_);

    const auto v = solve_thomas(sys);
    auto x = vector{};
    x.reserve(g.size());
    for (const auto i : std::views::iota(0ul, g.size())) x.push_back(g.x(i));

    return {x, v};
}

}  // namespace heat
