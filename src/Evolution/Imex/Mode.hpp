// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "Options/Options.hpp"

namespace imex {
/// IMEX implementations
enum class Mode {
  /// Use explicit stepping for IMEX terms.  The IMEX terms are still
  /// evolved separately from the explicit terms and will not be
  /// included in any time-stepper error estimation.
  Explicit,
  /// Solve the implicit equation using a nonlinear solver.
  Implicit,
  /// Solve a linearized version of the implicit equation.
  SemiImplicit,
};
}  // namespace imex

template <>
struct Options::create_from_yaml<imex::Mode> {
  template <typename Metavariables>
  static imex::Mode create(const Options::Option& options) {
    return create<void>(options);
  }
};

template <>
imex::Mode Options::create_from_yaml<imex::Mode>::create<void>(
    const Options::Option& options);
