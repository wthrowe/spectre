// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Time/TimeSteppers/Rk3PareschiSsp4.hpp"

namespace TimeSteppers {

size_t Rk3PareschiSsp4::order() const { return 3; }

size_t Rk3PareschiSsp4::error_estimate_order() const { return 2; }

double Rk3PareschiSsp4::stable_step() const { return 1.25637; }

size_t Rk3PareschiSsp4::imex_order() const { return 3; }

size_t Rk3PareschiSsp4::implicit_stage_order() const { return 0; }

namespace {
const double alpha = 0.24169426078821;
const double beta = 0.06042356519705;
const double eta = 0.12915286960590;
}  // namespace

const RungeKutta::ButcherTableau& Rk3PareschiSsp4::butcher_tableau() const {
  static const ButcherTableau tableau{
      // Substep times
      {0.0, 0.0, 1.0, 0.5},
      // Substep coefficients
      {{0.0},
       {0.0},
       {0.0, 0.0, 1.0},  // u1 = u0 + du0
       {0.0, 0.0, 1.0 / 4.0, 1.0 / 4.0}},
      // u2 = u0 + 1/4 du0 + 1/4 du1
      //    ~ 3/4 u0 + 1/4 (u1 + du1)  // c = 1
      // u3 = u0 + 1/6 du0 + 1/6 du1 + 2/3 du2
      //    ~ 1/3 u0 + 2/3 (u2 + du2)  // c = 1

      // Result coefficients
      {0.0, 0.0, 1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0},
      // Coefficients for the embedded method for generating an error measure.
      //
      // Not given in the reference.  Any set of coefficients of the form
      // (0, 0, x, x, 1 - 2 x) works.
      {0.0, 0.0, 0.0, 0.0, 1.0},
      // Dense output coefficient polynomials
      {{0.0, 1.0, -1.0},
       {},
       {0.0, 0.0, 1.0 / 6.0},
       {0.0, 0.0, 1.0 / 6.0},
       {0.0, 0.0, 2.0 / 3.0}}};
  return tableau;
}

const ImexRungeKutta::ImplicitButcherTableau&
Rk3PareschiSsp4::implicit_butcher_tableau() const {
  static const ImplicitButcherTableau tableau{
      {{0.0, alpha},  // This stage is zeroth-order
       {0.0, -alpha, alpha},
       {0.0, 0.0, 1 - alpha, alpha},
       {0.0, beta, eta, 0.5 - beta - eta - alpha, alpha}}};
  return tableau;
}
}  // namespace TimeSteppers

PUP::able::PUP_ID TimeSteppers::Rk3PareschiSsp4::my_PUP_ID = 0;  // NOLINT
