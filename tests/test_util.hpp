// SLAMCore-8G -- shared test scaffolding.
//
// A tiny assertion harness instead of a framework dependency: the project's
// whole premise is "git clone && make", so the test binary must build with
// nothing but a C++17 compiler.
#pragma once

#include <cmath>
#include <cstdio>
#include <string>

namespace sctest {

inline int& failures() {
  static int n = 0;
  return n;
}
inline int& checks() {
  static int n = 0;
  return n;
}

inline void report(const char* file, int line, const std::string& what, bool ok) {
  ++checks();
  if (ok) return;
  ++failures();
  std::printf("FAIL %s:%d  %s\n", file, line, what.c_str());
}

inline void near(const char* file, int line, const std::string& what, double got, double want,
                 double tol) {
  const bool ok = std::isfinite(got) && std::abs(got - want) <= tol;
  report(file, line, what + " (got " + std::to_string(got) + ", want " + std::to_string(want) +
                          " +-" + std::to_string(tol) + ")",
         ok);
}

inline void truth(const char* file, int line, const std::string& what, bool ok) {
  report(file, line, what, ok);
}

inline int finish(const char* name) {
  std::printf("%s: %d checks, %d failures\n", name, checks(), failures());
  return failures() == 0 ? 0 : 1;
}

}  // namespace sctest

#define SC_CHECK(cond) ::sctest::truth(__FILE__, __LINE__, #cond, (cond))
#define SC_NEAR(got, want, tol) ::sctest::near(__FILE__, __LINE__, #got, (got), (want), (tol))
#define SC_TEST_MAIN(name)            \
  int main() {                        \
    run();                            \
    return ::sctest::finish(name);    \
  }
