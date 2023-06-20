// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include "Framework/TestCreation.hpp"
#include "Framework/TestHelpers.hpp"
#include "Helpers/Time/TimeSteppers/RungeKutta.hpp"
#include "Helpers/Time/TimeSteppers/TimeStepperTestUtils.hpp"
#include "Time/TimeSteppers/Rk3PareschiSsp4.hpp"
#include "Time/TimeSteppers/TimeStepper.hpp"

SPECTRE_TEST_CASE("Unit.Time.TimeSteppers.Rk3PareschiSsp4", "[Unit][Time]") {
  const TimeSteppers::Rk3PareschiSsp4 stepper{};

  CHECK(stepper.order() == 3);
  CHECK(stepper.error_estimate_order() == 2);
  CHECK(stepper.number_of_substeps() == 5);
  CHECK(stepper.number_of_substeps_for_error() == 5);
  TestHelpers::RungeKutta::check_tableau(stepper);
  TestHelpers::RungeKutta::check_implicit_tableau(stepper, false);

  TimeStepperTestUtils::check_substep_properties(stepper);
  TimeStepperTestUtils::integrate_test(stepper, 3, 0, 1.0, 1.0e-9);
  TimeStepperTestUtils::integrate_test(stepper, 3, 0, -1.0, 1.0e-9);
  TimeStepperTestUtils::integrate_test_explicit_time_dependence(stepper, 3, 0,
                                                                -1.0, 1.0e-9);
  TimeStepperTestUtils::integrate_error_test(stepper, 3, 0, 1.0, 1.0e-8, 100,
                                             1.0e-4);
  TimeStepperTestUtils::integrate_error_test(stepper, 3, 0, -1.0, 1.0e-8, 100,
                                             1.0e-4);
  TimeStepperTestUtils::integrate_variable_test(stepper, 3, 0, 1.0e-9);
  TimeStepperTestUtils::stability_test(stepper);
  TimeStepperTestUtils::check_convergence_order(stepper, {10, 50});
  TimeStepperTestUtils::check_dense_output(stepper, 3_st);

  TimeStepperTestUtils::check_imex_convergence_order(stepper, {10, 50});

  TestHelpers::test_factory_creation<TimeStepper,
                                     TimeSteppers::Rk3PareschiSsp4>(
      "Rk3PareschiSsp4");
  test_serialization(stepper);
  test_serialization_via_base<TimeStepper, TimeSteppers::Rk3PareschiSsp4>();
  // test operator !=
  CHECK_FALSE(stepper != stepper);
}
