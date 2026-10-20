#include "app/app.hpp"
#include "core/errors.hpp"
#include "io/config_loader.hpp"
#include "io/writer.hpp"

#include <exception>
#include <print>

int main(int argc, char** argv) {
    if (argc != 2) {
        std::print(stderr, "Usage: {} <config.toml>\n", argv[0]);
        return 1;
    }

    try {
        const auto cfg = heat::load_config(argv[1]);
        const auto res = heat::solve_with_refinement(cfg);
        heat::write_csv(res, cfg.output_path_);

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
