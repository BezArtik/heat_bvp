// writer.hpp
#pragma once

#include <string>

namespace heat {

struct solve_result;
void write_csv(const solve_result& r, const std::string& path);

}  // namespace heat
