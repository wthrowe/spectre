// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include "Helpers/DataStructures/DataBox/TestHelpers.hpp"
#include "Evolution/Imex/Protocols/ImplicitSector.hpp"
#include "Evolution/Imex/Protocols/ImplicitSource.hpp"
#include "Evolution/Imex/Protocols/ImplicitSourceJacobian.hpp"
#include "Evolution/Imex/Tags/ImplicitHistory.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/Protocols/StaticReturnApplyable.hpp"

namespace {
struct Sector : tt::ConformsTo<imex::protocols::ImplicitSector> {
  using tensors = tmpl::list<>;

  struct source : tt::ConformsTo<imex::protocols::ImplicitSource>,
                  tt::ConformsTo<protocols::StaticReturnApplyable> {
    using return_tags = tmpl::list<>;
    using argument_tags = tmpl::list<>;
    static void apply();
  };

  struct source_jacobian
      : tt::ConformsTo<imex::protocols::ImplicitSourceJacobian>,
        tt::ConformsTo<protocols::StaticReturnApplyable> {
    using return_tags = tmpl::list<>;
    using argument_tags = tmpl::list<>;
    static void apply();
  };
};
}  // namespace

SPECTRE_TEST_CASE("Unit.Evolution.Imex.Tags.ImplicitHistory",
                  "[Unit][Evolution]") {
  TestHelpers::db::test_simple_tag<
      imex::Tags::ImplicitHistory<Sector>>("ImplicitHistory");
}
