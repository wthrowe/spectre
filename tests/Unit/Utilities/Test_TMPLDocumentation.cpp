// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cstddef>
#include <cstdint>

/// [include]
#include "Utilities/TMPL.hpp"
/// [include]

// We want the code to be nicely formatted for the documentation, not here.
// clang-format off
namespace {

template <typename T, typename U>
void assert_same() noexcept {
  static_assert(std::is_same_v<T, U>);
}

/// [example_declarations]
struct Type1;
struct Type2;
struct Type3;

// Usually, tmpl::list will be used in place of these.
template <typename...>
struct List1;
template <typename...>
struct List2;
template <typename...>
struct List3;
/// [example_declarations]

// Section: Containers

void integral_constant() noexcept {
/// [tmpl::integral_constant]
using T = tmpl::integral_constant<int, 3>;
static_assert(std::is_same_v<T::value_type, int>);
static_assert(std::is_same_v<T::type, T>);
static_assert(T::value == 3);

// At runtime only
CHECK(T{} == 3);
CHECK(T{}() == 3);
/// [tmpl::integral_constant]
//}

/// [tmpl::integral_constant::abbreviations]
static_assert(std::is_same_v<tmpl::int8_t<3>,
                             tmpl::integral_constant<int8_t, 3>>);
static_assert(std::is_same_v<tmpl::int16_t<3>,
                             tmpl::integral_constant<int16_t, 3>>);
assert_same<tmpl::int32_t<3>, tmpl::integral_constant<int32_t, 3>>();
assert_same<tmpl::int64_t<3>, tmpl::integral_constant<int64_t, 3>>();
}

static_assert(std::is_same_v<tmpl::uint8_t<3>,
                             tmpl::integral_constant<uint8_t, 3>>);
static_assert(std::is_same_v<tmpl::uint16_t<3>,
                             tmpl::integral_constant<uint16_t, 3>>);
static_assert(std::is_same_v<tmpl::uint32_t<3>,
                             tmpl::integral_constant<uint32_t, 3>>);
static_assert(std::is_same_v<tmpl::uint64_t<3>,
                             tmpl::integral_constant<uint64_t, 3>>);

static_assert(std::is_same_v<tmpl::size_t<3>,
                             tmpl::integral_constant<size_t, 3>>);
static_assert(std::is_same_v<tmpl::ptrdiff_t<3>,
                             tmpl::integral_constant<ptrdiff_t, 3>>);
static_assert(std::is_same_v<tmpl::bool_<true>,
                             tmpl::integral_constant<bool, true>>);
/// [tmpl::integral_constant::abbreviations]

/// [tmpl::list]
static_assert(not std::is_same_v<tmpl::list<Type1, Type2>,
                                 tmpl::list<Type2, Type1>>);
/// [tmpl::list]

/// [tmpl::type_]
static_assert(std::is_same_v<tmpl::type_<Type1>::type, Type1>);
/// [tmpl::type_]

// Section: Constants

/// [tmpl::empty_sequence]
static_assert(std::is_same_v<tmpl::empty_sequence, tmpl::list<>>);
/// [tmpl::empty_sequence]

/// [tmpl::false_type]
static_assert(std::is_same_v<tmpl::false_type, tmpl::bool_<false>>);
/// [tmpl::false_type]

/// [tmpl::no_such_type_]
static_assert(std::is_same_v<tmpl::index_of<List1<>, Type1>,
                             tmpl::no_such_type_>);
/// [tmpl::no_such_type_]

/// [tmpl::true_type]
static_assert(std::is_same_v<tmpl::true_type, tmpl::bool_<true>>);
/// [tmpl::true_type]

// Constructor-like functions for lists

/// [tmpl::filled_list]
static_assert(std::is_same_v<tmpl::filled_list<Type1, 3, List1>,
                             List1<Type1, Type1, Type1>>);
static_assert(std::is_same_v<tmpl::filled_list<Type1, 3>,
                             tmpl::list<Type1, Type1, Type1>>);
/// [tmpl::filled_list]

/// [tmpl::integral_list]
static_assert(std::is_same_v<tmpl::integral_list<int, 3, 2, 1>,
                             tmpl::list<tmpl::integral_constant<int, 3>,
                                        tmpl::integral_constant<int, 2>,
                                        tmpl::integral_constant<int, 1>>>);
/// [tmpl::integral_list]

// Section: Functions for querying lists

/// [tmpl::at_c]
static_assert(std::is_same_v<tmpl::at_c<List1<Type1, Type2, Type3>, 0>, Type1>);
/// [tmpl::at_c]

/// [tmpl::back]
static_assert(std::is_same_v<tmpl::back<List1<Type1, Type2, Type3>>, Type3>);
/// [tmpl::back]

/// [tmpl::front]
static_assert(std::is_same_v<tmpl::front<List1<Type1, Type2, Type3>>, Type1>);
/// [tmpl::front]

/// [tmpl::index_of]
static_assert(std::is_same_v<tmpl::index_of<List1<Type1, Type2, Type3>, Type3>,
                             tmpl::size_t<2>>);
static_assert(std::is_same_v<tmpl::index_of<List1<Type1, Type3, Type3>, Type3>,
                             tmpl::size_t<1>>);
static_assert(std::is_same_v<tmpl::index_of<List1<Type1>, Type2>,
                             tmpl::no_such_type_>);
/// [tmpl::index_of]

/// [tmpl::list_contains]
static_assert(std::is_same_v<tmpl::list_contains<List1<Type1, Type2>, Type1>,
                             tmpl::integral_constant<bool, true>>);
static_assert(tmpl::list_contains_v<List1<Type1, Type2>, Type1>);
static_assert(not tmpl::list_contains_v<List1<Type2, Type2>, Type1>);
/// [tmpl::list_contains]

/// [tmpl::size]
static_assert(std::is_same_v<tmpl::size<List1<Type1, Type1>>,
                             tmpl::integral_constant<unsigned int, 2>>);
/// [tmpl::size]

// Section: Functions producing lists from other lists

/// [tmpl::append]
static_assert(
  std::is_same_v<tmpl::append<List1<Type1, Type2>, List2<>, List2<Type2>>,
                 List1<Type1, Type2, Type2>>);
/// [tmpl::append]

/// [tmpl::append::bug]
static_assert(
  std::is_same_v<tmpl::append<List1<Type1, Type2>, List2<>, Type2>,
                 tmpl::list<>>);
/// [tmpl::append::bug]

/// [tmpl::clear]
static_assert(std::is_same_v<tmpl::clear<List1<Type1>>, List1<>>);
/// [tmpl::clear]

/// [tmpl::erase_c]
static_assert(std::is_same_v<tmpl::erase_c<List1<Type1, Type2, Type3>, 1>,
                             List1<Type1, Type3>>);
/// [tmpl::erase_c]

/// [tmpl::join]
static_assert(
  std::is_same_v<tmpl::join<List3<List1<Type1, Type2>, List2<>, List2<Type2>>>,
                 List1<Type1, Type2, Type2>>);
/// [tmpl::join]

/// [tmpl::join::bug]
static_assert(
  std::is_same_v<tmpl::join<List3<List1<Type1, Type2>, List2<>, Type2>>,
                 tmpl::list<>>);
/// [tmpl::join::bug]

/// [tmpl::pop_back]
static_assert(std::is_same_v<tmpl::pop_back<List1<Type1, Type2, Type3>>,
                             List1<Type1, Type2>>);
static_assert(
  std::is_same_v<tmpl::pop_back<List1<Type1, Type2, Type3>,
                                tmpl::integral_constant<unsigned int, 2>>,
                 List1<Type1>>);
/// [tmpl::pop_back]

/// [tmpl::pop_front]
static_assert(std::is_same_v<tmpl::pop_front<List1<Type1, Type2, Type3>>,
                             List1<Type2, Type3>>);
static_assert(
  std::is_same_v<tmpl::pop_front<List1<Type1, Type2, Type3>,
                                 tmpl::integral_constant<unsigned int, 2>>,
                 List1<Type3>>);
/// [tmpl::pop_front]

/// [tmpl::push_back]
static_assert(std::is_same_v<tmpl::push_back<List1<Type1>, Type2, Type3>,
                             List1<Type1, Type2, Type3>>);
/// [tmpl::push_back]

/// [tmpl::push_front]
static_assert(std::is_same_v<tmpl::push_front<List1<Type1>, Type2, Type3>,
                             List1<Type2, Type3, Type1>>);
/// [tmpl::push_front]

/// [tmpl::split_at]
static_assert(
  std::is_same_v<tmpl::split_at<List1<Type1, Type2, Type3>,
                                tmpl::integral_constant<unsigned int, 2>>,
                 List1<List1<Type1, Type2>, List1<Type3>>>);
/// [tmpl::split_at]

// Mathematical functions

/// [math_bitwise]
static_assert(
  std::is_same_v<tmpl::complement<tmpl::uint8_t<0b10001111>>::type,
                                  tmpl::uint8_t<0b01110000>>);
static_assert(
  std::is_same_v<tmpl::bitand_<tmpl::uint8_t<0b00111011>,
                               tmpl::uint8_t<0b01010110>>::type,
                               tmpl::uint8_t<0b00010010>>);
static_assert(
  std::is_same_v<tmpl::bitor_<tmpl::uint8_t<0b01100011>,
                              tmpl::uint8_t<0b10100111>>::type,
                              tmpl::uint8_t<0b11100111>>);
static_assert(
  std::is_same_v<tmpl::bitxor_<tmpl::uint8_t<0b11000011>,
                               tmpl::uint8_t<0b00000110>>::type,
                               tmpl::uint8_t<0b11000101>>);
/// [math_bitwise]

} // namespace
// clang-format on

SPECTRE_TEST_CASE("Unit.Utilities.TMPL.Documentation", "[Unit][Utilities]") {
  integral_constant();
}
