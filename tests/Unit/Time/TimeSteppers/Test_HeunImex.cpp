// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include "Framework/TestCreation.hpp"
#include "Framework/TestHelpers.hpp"
#include "Helpers/Time/TimeSteppers/TimeStepperTestUtils.hpp"
#include "Time/TimeSteppers/HeunImex.hpp"
#include "Time/TimeSteppers/TimeStepper.hpp"

SPECTRE_TEST_CASE("Unit.Time.TimeSteppers.HeunImex", "[Unit][Time]") {
  const TimeSteppers::HeunImex stepper{};
  TimeStepperTestUtils::check_substep_properties(stepper);
  TimeStepperTestUtils::integrate_test(stepper, 2, 0, 1., 1e-6);
  TimeStepperTestUtils::integrate_test(stepper, 2, 0, -1., 1e-6);
  TimeStepperTestUtils::integrate_test_explicit_time_dependence(stepper, 2, 0,
                                                                -1.0, 1.0e-6);
  TimeStepperTestUtils::integrate_error_test(stepper, 2, 0, 1.0, 1.0e-5, 100,
                                             1.0e-3);
  TimeStepperTestUtils::integrate_error_test(stepper, 2, 0, -1.0, 1.0e-5, 100,
                                             1.0e-3);
  TimeStepperTestUtils::integrate_variable_test(stepper, 2, 0, 1e-6);
  TimeStepperTestUtils::stability_test(stepper);
  TimeStepperTestUtils::check_convergence_order(stepper);
  TimeStepperTestUtils::check_dense_output(stepper, 2);

  CHECK(stepper.order() == 2);
  CHECK(stepper.error_estimate_order() == 1);

  TestHelpers::test_factory_creation<TimeStepper, TimeSteppers::HeunImex>(
      "HeunImex");
  test_serialization(stepper);
  test_serialization_via_base<TimeStepper, TimeSteppers::HeunImex>();
  // test operator !=
  CHECK_FALSE(stepper != stepper);
}
