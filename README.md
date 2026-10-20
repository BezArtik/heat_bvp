# heat_bvp

A C++23 command-line utility for solving 1D boundary value problems for the
stationary heat equation with piecewise coefficients, using the finite volume
(balance) method and the Thomas algorithm.

--------------------------------------------------------------------------------
1. MATHEMATICAL FORMULATION
--------------------------------------------------------------------------------

The solver targets the following boundary value problem:

    d/dx ( k(x) * du/dx ) - q(x) * u(x) = -f(x),    x in (0, 1)
    u(0) = mu1,   u(1) = mu2

The coefficients k(x), q(x), and f(x) are piecewise functions with a
discontinuity at a given point xi.
The solution u(x) and the heat flux w(x) = -k(x)*u'(x) are continuous.

The problem is discretized using a conservative finite-difference scheme
(balance method) on a uniform grid, resulting in a tridiagonal linear system
which is solved using the O(n) Thomas algorithm.

--------------------------------------------------------------------------------
2. PROJECT STRUCTURE
--------------------------------------------------------------------------------

    heat_bvp/
    |-- CMakeLists.txt          # Root build file
    |-- CMakePresets.json      # Compiler presets (clang-release, clang-debug)
    |
    |-- include/                # Header files
    |   |-- app/                # Orchestration and refinement logic
    |   |-- core/               # Types, errors, config structures
    |   |-- io/                 # Config loader, CSV writer
    |   |-- math/               # Grid, quadrature, tridiagonal solver, scheme, BVP solver
    |
    |-- src/                    # Source files
        |-- app/
        |-- io/
        |-- problems/
        |-- main.cpp            # Entry point

--------------------------------------------------------------------------------
3. REQUIREMENTS
--------------------------------------------------------------------------------

C++ Compiler: Clang 16+, GCC 13+, or MSVC 19.30+
CMake: 3.24+

Dependencies are fetched automatically via CMake FetchContent:
- toml++       (TOML parsing)
- ExprTk       (Analytical expressions)

--------------------------------------------------------------------------------
4. BUILDING
--------------------------------------------------------------------------------

Configure and build the project in release mode:

    cmake --preset clang-release
    cmake --build --preset clang-release

The executable will be located at:

    out/build/clang-release/heat_app

For a debug build with sanitizers:

    cmake --preset clang-debug
    cmake --build --preset clang-debug

--------------------------------------------------------------------------------
5. USAGE
--------------------------------------------------------------------------------

Run the solver by passing a TOML configuration file:

    ./out/build/clang-release/heat_app configs/test.toml

The program automatically refines the grid until the requested accuracy (eps)
is achieved and prints a summary to stderr.

    Example output:
    n = 640
    iterations = 7
    eps1 = 3.644459216634033e-07
    max at x = 0.578125

The detailed numerical results (node index, x, numerical solution,
reference solution, and error) are saved to the path specified in the config.

--------------------------------------------------------------------------------
6. CONFIGURATION FORMAT (.toml)
--------------------------------------------------------------------------------

    [problem]
    a  = 0.0
    b  = 1.0
    xi = 0.5773502691896257
    mu1 = 2.0
    mu2 = 1.0
    boundary_type = "dirichlet"
    // Optional analytical expression for test mode (leave empty for main mode)
    analytical = "(x <= xi) ? (...) : (...)"

    [problem.left]
    k = "exp(xi^2)"
    q = "xi^2"
    f = "xi^2 - 1.0"

    [problem.right]
    k = "xi^2"
    q = "1.0 + xi^4"
    f = "1.0"

    [solver]
    n     = 10
    eps   = 5e-7
    max_n = 100000

    [output]
    path = "labwork/output/test.csv"
    print_every = 1

--------------------------------------------------------------------------------
8. LICENSE
--------------------------------------------------------------------------------

MIT
