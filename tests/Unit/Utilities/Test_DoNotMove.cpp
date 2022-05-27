// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include "Utilities/DoNotMove.hpp"

namespace {
enum class Source { Default, Copied, Moved };
struct RecordAssignment {
  RecordAssignment() = default;
  RecordAssignment(const RecordAssignment&) = delete;
  RecordAssignment(RecordAssignment&&) = delete;

  RecordAssignment& operator=(const RecordAssignment&) {
    source = Source::Copied;
    return *this;
  }

  RecordAssignment& operator=(RecordAssignment&&) {
    source = Source::Moved;
    return *this;
  }

  Source source = Source::Default;
};

RecordAssignment make_rvalue() { return RecordAssignment(); }
RecordAssignment& make_lvalue() {
  static RecordAssignment result;
  return result;
}
}  // namespace

SPECTRE_TEST_CASE("Unit.Utilities.do_not_move", "[Unit][Utilities]") {
  RecordAssignment rvalue_default{};
  rvalue_default = make_rvalue();
  RecordAssignment rvalue_dnm{};
  rvalue_dnm = do_not_move(make_rvalue());
  RecordAssignment lvalue_default{};
  lvalue_default = make_lvalue();
  RecordAssignment lvalue_dnm{};
  lvalue_dnm = do_not_move(make_lvalue());
  CHECK(rvalue_default.source == Source::Moved);
  CHECK(rvalue_dnm.source == Source::Copied);
  CHECK(lvalue_default.source == Source::Copied);
  CHECK(lvalue_dnm.source == Source::Copied);
}
