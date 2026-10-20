// writer.cpp
#include "io/writer.hpp"

#include "app/app.hpp"
#include "core/errors.hpp"

#include <charconv>
#include <format>
#include <fstream>
#include <string>

namespace heat {

namespace {

constexpr std::size_t flush_thresh = 64 * 1024;

void append_num(std::string& b, double v) {
    std::array<char, 32> tmp;
    auto [p, _] = std::to_chars(tmp.begin(), tmp.end(), v);
    b.append(tmp.begin(), p);
}

void append_int(std::string& b, std::size_t v) {
    std::array<char, 24> tmp;
    auto [p, _] = std::to_chars(tmp.begin(), tmp.end(), v);
    b.append(tmp.begin(), p);
}

}  // namespace

void write_csv(const solve_result& r, const std::string& path) {
    auto file = std::ofstream{path};
    if (!file.is_open()) throw solver_error{std::format("cannot open '{}'", path)};

    auto buf = std::string{};
    buf.reserve(flush_thresh + 256);
    buf += r.is_test_ ? "i,x,u,v,err\n" : "i,x,v,v2,err\n";

    for (std::size_t i = 0; i < r.x_.size(); ++i) {
        append_int(buf, i);
        buf += ',';
        append_num(buf, r.x_[i]);
        buf += ',';
        append_num(buf, r.v_[i]);
        buf += ',';
        append_num(buf, r.reference_[i]);
        buf += ',';
        append_num(buf, r.error_[i]);
        buf += '\n';

        if (buf.size() >= flush_thresh) {
            file.write(buf.data(), buf.size());
            buf.clear();
        }
    }
    if (!buf.empty()) file.write(buf.data(), buf.size());
}
}  // namespace heat
