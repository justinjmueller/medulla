/**
 * @file test_systematic_clip.cc
 * @brief Standalone check of sys::Systematic::clip and the clip-range
 * validation in the Systematic constructor.
 * @details std::clamp returns a NaN argument unchanged, so a NaN weight used
 * to pass straight through the clip and, via TH2::Fill, permanently turn the
 * affected histogram bin into NaN. The policy is that a NaN weight is always
 * replaced by 1.0 (no reweighting), whatever the configured range. The clip
 * range must be empty (no clipping) or exactly two values [min, max]; any other
 * length used to be silently accepted (one value disabled clipping, more than
 * two ignored the extras) and is now rejected with std::invalid_argument.
 *
 * Build (with the environment from setup.sh), from the build directory layout
 * used by the top-level CMake project:
 *   B=../../build
 *   g++ -std=c++17 -Wall -O2 \
 *       -I../include -I../../shared/include -I$B/_deps/tomlplusplus-src/include \
 *       -I$(root-config --incdir) \
 *       test_systematic_clip.cc ../src/systematic.cc \
 *       -o test_systematic_clip \
 *       -L$B/shared -lshared -Wl,-rpath,$PWD/$B/shared $(root-config --libs)
 *
 * Run:
 *   ./test_systematic_clip
 *
 * Exit status is 0 if every check passes and 1 otherwise.
 */
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include <unistd.h>

#include "configuration.h"
#include "systematic.h"

static int failures = 0;
#define CHECK(cond, msg) do { if(cond) std::printf("ok:   %s\n", msg); else { ++failures; std::printf("FAIL: %s\n", msg); } } while(0)

// Holds a parsed configuration and its [[sys]] tables. Heap-allocated and
// never moved, because the tables point into the parent configuration.
struct Case
{
    cfg::ConfigurationTable config;
    std::vector<cfg::ConfigurationTable> subs;
};

static Case * make_case(const std::string & clip_line)
{
    std::string path = "/tmp/test_systematic_clip_" + std::to_string(getpid()) + ".toml";
    {
        std::ofstream out(path);
        out << "[general]\nname = \"clip_test\"\n\n[[sys]]\nname = \"t\"\ntype = \"multisim\"\nindex = 0\n" << clip_line << "\n";
    }
    Case * c = new Case();
    c->config.set_config(path);
    c->subs = c->config.get_subtables("sys");
    std::remove(path.c_str());
    return c;
}

static const double nan_value = std::numeric_limits<double>::quiet_NaN();
static const double inf_value = std::numeric_limits<double>::infinity();

static bool throws_invalid(const std::string & clip_line, std::vector<double> default_clip)
{
    Case * c = make_case(clip_line);
    try { sys::Systematic s(c->subs[0], nullptr, default_clip); }
    catch(const std::invalid_argument &) { return true; }
    return false;
}

int main()
{
    // Range from the systematic own "clip" field.
    {
        Case * c = make_case("clip = [0.0, 10.0]");
        sys::Systematic s(c->subs[0], nullptr, {});
        CHECK(s.clip(nan_value) == 1.0, "NaN -> 1.0 with clip [0, 10]");
        CHECK(s.clip(-nan_value) == 1.0, "negative-signed NaN -> 1.0");
        CHECK(s.clip(std::nan("")) == 1.0, "std::nan NaN -> 1.0");
        CHECK(s.clip(inf_value) == 10.0, "+inf still clipped to max");
        CHECK(s.clip(-inf_value) == 0.0, "-inf still clipped to min");
        CHECK(s.clip(-1.0) == 0.0, "below range clipped to min");
        CHECK(s.clip(42.0) == 10.0, "above range clipped to max");
        CHECK(s.clip(5.0) == 5.0, "in-range weight unchanged");
    }

    // NaN is 1.0 even when the configured range excludes 1.0.
    {
        Case * c = make_case("clip = [2.0, 3.0]");
        sys::Systematic s(c->subs[0], nullptr, {});
        CHECK(s.clip(nan_value) == 1.0, "NaN -> 1.0 even if range [2, 3] excludes 1");
        CHECK(s.clip(1.5) == 2.0, "ordinary weight still clipped to [2, 3]");
    }

    // No clip configured anywhere: no clipping, but NaN is still replaced.
    {
        Case * c = make_case("");
        sys::Systematic s(c->subs[0], nullptr, {});
        CHECK(s.clip(nan_value) == 1.0, "NaN -> 1.0 with no clip configured");
        CHECK(s.clip(42.0) == 42.0 && s.clip(-1.0) == -1.0, "finite weights unchanged with no clip");
        CHECK(std::isinf(s.clip(inf_value)), "inf passes through with no clip");
    }

    // General default applies unless the systematic overrides it.
    {
        Case * c = make_case("");
        sys::Systematic s(c->subs[0], nullptr, {0.0, 10.0});
        CHECK(s.clip(nan_value) == 1.0, "NaN -> 1.0 with default clip");
        CHECK(s.clip(42.0) == 10.0, "default clip range applied");
        Case * d = make_case("clip = [0.0, 2.0]");
        sys::Systematic t(d->subs[0], nullptr, {0.0, 10.0});
        CHECK(t.clip(5.0) == 2.0, "own clip overrides default");
    }

    // An explicitly empty clip array means no clipping.
    {
        Case * c = make_case("clip = []");
        bool ok = true;
        try { sys::Systematic s(c->subs[0], nullptr, {0.0, 10.0}); CHECK(s.clip(42.0) == 42.0, "empty own clip array disables clipping"); }
        catch(...) { ok = false; }
        CHECK(ok, "empty own clip array is accepted");
    }

    // Length validation.
    CHECK(throws_invalid("clip = [1.0]", {}), "own clip with 1 value throws");
    CHECK(throws_invalid("clip = [0.0, 1.0, 2.0]", {}), "own clip with 3 values throws");
    CHECK(throws_invalid("", {5.0}), "default clip with 1 value throws");
    CHECK(throws_invalid("", {0.0, 1.0, 2.0}), "default clip with 3 values throws");
    CHECK(!throws_invalid("clip = [0.0, 1.0]", {}), "own clip with 2 values does not throw");

    std::printf("%s (%d failure%s)\n", failures ? "FAILED" : "PASSED", failures, failures == 1 ? "" : "s");
    return failures ? 1 : 0;
}
