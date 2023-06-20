// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "Options/String.hpp"
#include "Time/TimeSteppers/ImexRungeKutta.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/TMPL.hpp"

namespace TimeSteppers {
/*!
 * \ingroup TimeSteppersGroup
 * \brief A fourth-order strong-stability-preserving Runge-Kutta
 * method with IMEX support.
 *
 * The method as published has four stages, but is implemented with
 * five as a way to convert it to an EDIRK method.
 *
 * The coefficients are given as IMEX-SSP3(4,3,3) in \cite Pareschi2005.
 *
 * Using this time stepper in a non-IMEX simulation is not
 * recommended, as it performs two unused RHS evaluations.  When using
 * IMEX it performs one extra evaluation because there are more
 * implicit steps than explicit.
 *
 * The implicit portion is L-stable.
 *
 * The CFL factor/stable step size is 1.25637.
 */
class Rk3PareschiSsp4 : public ImexRungeKutta {
 public:
  using options = tmpl::list<>;
  static constexpr Options::String help = {
      "A 3rd-order 4 stage Runge-Kutta scheme devised by Pareschi and Russo."};

  Rk3PareschiSsp4() = default;
  Rk3PareschiSsp4(const Rk3PareschiSsp4&) = default;
  Rk3PareschiSsp4& operator=(const Rk3PareschiSsp4&) = default;
  Rk3PareschiSsp4(Rk3PareschiSsp4&&) = default;
  Rk3PareschiSsp4& operator=(Rk3PareschiSsp4&&) = default;
  ~Rk3PareschiSsp4() override = default;

  size_t order() const override;

  size_t error_estimate_order() const override;

  double stable_step() const override;

  size_t imex_order() const override;

  size_t implicit_stage_order() const override;

  WRAPPED_PUPable_decl_template(Rk3PareschiSsp4);  // NOLINT

  explicit Rk3PareschiSsp4(CkMigrateMessage* /*unused*/) {}

  const ButcherTableau& butcher_tableau() const override;

  const ImplicitButcherTableau& implicit_butcher_tableau() const override;
};

inline bool constexpr operator==(const Rk3PareschiSsp4& /*lhs*/,
                                 const Rk3PareschiSsp4& /*rhs*/) {
  return true;
}

inline bool constexpr operator!=(const Rk3PareschiSsp4& lhs,
                                 const Rk3PareschiSsp4& rhs) {
  return not(lhs == rhs);
}
}  // namespace TimeSteppers
