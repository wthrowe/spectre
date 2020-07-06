// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <type_traits>

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

#define HAS_LAZY_VERSION(name) \
  assert_same<tmpl::wrap<lazy_test_arguments, tmpl::lazy::name>::type, \
              tmpl::wrap<lazy_test_arguments, tmpl::name>>()

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

template <typename>
struct Wrapper;

template <typename... T>
struct make_list1 {
  using type = List1<T...>;
};
/// [example_declarations]

namespace metafunctions {
/// [metafunctions:lazy]
template <typename T>
struct lazy_add_list {
  using type = List1<T, int>;
};
/// [metafunctions:lazy]

/// [metafunctions:eager]
template <typename T>
using eager_add_list = List1<T, int>;
/// [metafunctions:eager]

/// [metafunctions:eager_from_lazy]
template <typename T>
using eager_add_list2 = typename lazy_add_list<T>::type;
/// [metafunctions:eager_from_lazy]

/// [metafunctions:lazy_call_metafunction]
template <typename Func, typename Arg>
struct lazy_call_func_impl : Func {};
template <template <typename> typename Func, typename DummyArg, typename Arg>
struct lazy_call_func_impl<Func<DummyArg>, Arg> : Func<Arg> {};

template <typename Func>
using lazy_call_func = List1<typename lazy_call_func_impl<Func, double>::type,
                             typename lazy_call_func_impl<Func, char>::type>;
/// [metafunctions:lazy_call_metafunction]

/// [metafunctions:eager_call_metafunction]
template <typename Func, typename Arg>
struct eager_call_func_impl {
  using type = Func;
};
template <template <typename> typename Func, typename DummyArg, typename Arg>
struct eager_call_func_impl<Func<DummyArg>, Arg> {
  using type = Func<Arg>;
};

template <typename Func>
using eager_call_func = List1<typename eager_call_func_impl<Func, double>::type,
                              typename eager_call_func_impl<Func, char>::type>;
/// [metafunctions:eager_call_metafunction]

/// [metafunctions:evens]
template <typename L>
using evens = tmpl::filter<
  L, tmpl::equal_to<tmpl::modulo<tmpl::_1, tmpl::integral_constant<int, 2>>,
                    tmpl::integral_constant<int, 0>>>;
/// [metafunctions:evens]

/// [metafunctions:maybe_first]
template <typename L>
using maybe_first = tmpl::apply<tmpl::apply<
  tmpl::if_<tmpl::size<L>,
            tmpl::defer<tmpl::bind<tmpl::front, tmpl::pin<L>>>,
            tmpl::no_such_type_>>>;
/// [metafunctions:maybe_first]

/// [metafunctions:factorial]
template <typename N>
using factorial =
  tmpl::fold<tmpl::range<typename N::value_type, 1, N::value + 1>,
             tmpl::integral_constant<typename N::value_type, 1>,
             tmpl::times<tmpl::_state, tmpl::_element>>;
/// [metafunctions:factorial]

/// [metafunctions:multiplication_table]
template <typename N>
using multiplication_table =
  tmpl::transform<
    tmpl::range<typename N::value_type, 1, N::value + 1>,
    tmpl::lazy::transform<
      tmpl::pin<tmpl::range<typename N::value_type, 1, N::value + 1>>,
      tmpl::defer<tmpl::times<tmpl::_1, tmpl::parent<tmpl::_1>>>>>;
/// [metafunctions:multiplication_table]

/// [metafunctions:column_with_zeros]
template <typename Lists, typename Column>
using column_with_zeros =
  tmpl::transform<
    Lists,
    tmpl::bind<
      tmpl::apply,
      tmpl::if_<
        tmpl::greater<tmpl::bind<tmpl::size, tmpl::_1>, Column>,
        tmpl::defer<  // avoid out-of-range call to `at`
          tmpl::always<
            tmpl::parent<
              tmpl::bind<tmpl::at, tmpl::_1, Column>>>>,
        tmpl::size_t<0>>>>;
/// [metafunctions:column_with_zeros]

/// [metafunctions:factorial_recursion]
template <typename N>
using factorial_recursion =
  tmpl::apply<
    tmpl::bind<tmpl::apply, tmpl::_1,
               tmpl::_1, tmpl::integral_constant<typename N::value_type, 1>, N>,
    tmpl::bind<  // recursive metalambda starts here
      tmpl::apply,
      tmpl::if_<
        tmpl::_3,
        tmpl::defer<  // prevent speculative recursion
          tmpl::always<
            tmpl::parent<
              tmpl::bind<tmpl::apply, tmpl::_1,
                         tmpl::_1,
                         tmpl::times<tmpl::_2, tmpl::_3>,
                         tmpl::prev<tmpl::_3>>>>>,
        tmpl::_2>>>;
/// [metafunctions:factorial_recursion]

/// [metafunctions:primes]
template <typename N>
using zero = tmpl::integral_constant<typename N::value_type, 0>;

template <typename Sequence, typename N, typename NewType>
using replace_at =
  tmpl::append<tmpl::front<tmpl::split_at<Sequence, N>>,
               tmpl::list<NewType>,
               tmpl::pop_front<tmpl::back<tmpl::split_at<Sequence, N>>>>;

template <typename Start, typename End>
using range_from_types =
  tmpl::range<typename Start::value_type, Start::value, End::value>;

template <typename N>
using primes =
  tmpl::filter<
    tmpl::fold<
      tmpl::range<typename N::value_type, 2, N::value>,
      tmpl::push_front<
        tmpl::range<typename N::value_type, 2, N::value>, zero<N>, zero<N>>,
      tmpl::bind<
        tmpl::apply,
        tmpl::if_<
          tmpl::or_<
            std::is_same<
              tmpl::bind<tmpl::at, tmpl::_state, tmpl::_element>, zero<N>>,
            tmpl::greater<tmpl::times<tmpl::_element, tmpl::_element>, N>>,
          tmpl::defer<  // needed to match other branch (don't apply<_state>)
            tmpl::always<tmpl::parent<tmpl::_state>>>,
          tmpl::defer<  // prevent math on argument placeholders
            tmpl::always<
              tmpl::parent<
                tmpl::lazy::fold<
                  tmpl::bind<
                    range_from_types,
                    tmpl::_element,
                    tmpl::next<tmpl::divides<tmpl::prev<N>, tmpl::_element>>>,
                  tmpl::_state,
                  tmpl::defer<  // passed as a closure to fold
                    tmpl::bind<
                      replace_at,
                      tmpl::_state,
                      tmpl::times<tmpl::_element, tmpl::parent<tmpl::_element>>,
                      zero<N>>>>>>>>>>,
    tmpl::not_equal_to<tmpl::_1, zero<N>>>;
/// [metafunctions:primes]
}  // namespace metafunctions

SPECTRE_TEST_CASE("Unit.Utilities.TMPL.Documentation", "[Unit][Utilities]") {
{
using namespace metafunctions;

/// [metafunctions:agreement]
assert_same<eager_add_list<double>, eager_add_list2<double>>();
/// [metafunctions:agreement]

/// [metafunctions:lazy_call_metafunction_assert]
struct Dummy;

assert_same<lazy_call_func<lazy_add_list<Dummy>>,
            List1<List1<double, int>, List1<char, int>>>();
/// [metafunctions:lazy_call_metafunction_assert]

/// [metafunctions:eager_call_metafunction_assert]
assert_same<eager_call_func<eager_add_list<Dummy>>,
            List1<List1<Dummy, int>, List1<Dummy, int>>>();
/// [metafunctions:eager_call_metafunction_assert]

/// [metafunctions:evens:asserts]
assert_same<evens<tmpl::integral_list<int, 1, 1, 2, 3, 5, 8, 13>>,
            tmpl::integral_list<int, 2, 8>>();
/// [metafunctions:evens:asserts]

/// [metafunctions:maybe_first:asserts]
assert_same<maybe_first<tmpl::list<Type1>>, Type1>();
assert_same<maybe_first<tmpl::list<Type1, Type2>>, Type1>();
assert_same<maybe_first<tmpl::list<>>, tmpl::no_such_type_>();
/// [metafunctions:maybe_first:asserts]

/// [metafunctions:factorial:asserts]
assert_same<factorial<tmpl::size_t<5>>, tmpl::size_t<120>>();
/// [metafunctions:factorial:asserts]

/// [metafunctions:multiplication_table:asserts]
assert_same<multiplication_table<tmpl::size_t<5>>,
            tmpl::list<tmpl::integral_list<size_t, 1, 2, 3, 4, 5>,
                       tmpl::integral_list<size_t, 2, 4, 6, 8, 10>,
                       tmpl::integral_list<size_t, 3, 6, 9, 12, 15>,
                       tmpl::integral_list<size_t, 4, 8, 12, 16, 20>,
                       tmpl::integral_list<size_t, 5, 10, 15, 20, 25>>>();
/// [metafunctions:multiplication_table:asserts]

/// [metafunctions:factorial_recursion:asserts]
assert_same<factorial_recursion<tmpl::size_t<5>>, tmpl::size_t<120>>();
/// [metafunctions:factorial_recursion:asserts]

/// [metafunctions:column_with_zeros:asserts]
assert_same<
  column_with_zeros<
    tmpl::list<tmpl::integral_list<size_t, 11, 12, 13>,
               tmpl::integral_list<size_t, 21, 22, 23, 24, 25>,
               tmpl::integral_list<size_t, 31, 32, 33, 34>>,
    tmpl::size_t<3>>,
  tmpl::integral_list<size_t, 0, 24, 34>>();
/// [metafunctions:column_with_zeros:asserts]

/// [metafunctions:primes:asserts]
assert_same<
  primes<tmpl::size_t<100>>,
  tmpl::integral_list<size_t, 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41,
                      43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97>>();
/// [metafunctions:primes:asserts]
}

/// [tmpl::args]
static_assert(not std::is_same_v<tmpl::_1, tmpl::args<0>>);
static_assert(not std::is_same_v<tmpl::_2, tmpl::args<1>>);
static_assert(std::is_same_v<tmpl::_3, tmpl::args<2>>);
static_assert(std::is_same_v<tmpl::_4, tmpl::args<3>>);
/// [tmpl::args]

/// [tmpl::args:eval]
assert_same<tmpl::apply<tmpl::_1, Type1, Type2>, Type1>();
assert_same<tmpl::apply<tmpl::_2, Type1, Type2>, Type2>();
assert_same<tmpl::apply<tmpl::args<0>, Type1, Type2>, Type1>();
/// [tmpl::args:eval]

/// [metalambda_lazy]
assert_same<tmpl::apply<make_list1<tmpl::_1, tmpl::_2>, Type1, Type2>,
            List1<Type1, Type2>>();
/// [metalambda_lazy]

/// [tmpl::bind]
assert_same<tmpl::apply<tmpl::bind<List1, tmpl::_1, tmpl::_2>, Type1, Type2>,
            List1<Type1, Type2>>();
/// [tmpl::bind]

/// [tmpl::pin]
assert_same<tmpl::apply<tmpl::pin<List1<Type1>>>, List1<Type1>>();
// Error: List1 is not a lazy metafunction
// assert_same<tmpl::apply<List1<Type1>>, List1<Type1>>();
assert_same<tmpl::apply<tmpl::pin<tmpl::_1>, Type1>, tmpl::_1>();
/// [tmpl::pin]

/// [tmpl::defer]
assert_same<
  tmpl::apply<tmpl::apply<tmpl::defer<tmpl::always<tmpl::_1>>,
                          Type1>,
              Type2>,
  Type2>();
/// [tmpl::defer]

/// [tmpl::parent]
assert_same<
  tmpl::apply<tmpl::apply<tmpl::defer<tmpl::always<tmpl::parent<tmpl::_1>>>,
                          Type1>,
              Type2>,
  Type1>();
/// [tmpl::parent]

/// [metalambda_constant]
assert_same<tmpl::apply<Type1>, Type1>();
/// [metalambda_constant]

// Section: Containers

/// [tmpl::integral_constant]
using T = tmpl::integral_constant<int, 3>;
assert_same<T::value_type, int>();
assert_same<T::type, T>();
static_assert(T::value == 3);

// At runtime only
CHECK(T{} == 3);
CHECK(T{}() == 3);
/// [tmpl::integral_constant]

/// [tmpl::integral_constant::abbreviations]
assert_same<tmpl::int8_t<3>, tmpl::integral_constant<int8_t, 3>>();
assert_same<tmpl::int16_t<3>, tmpl::integral_constant<int16_t, 3>>();
assert_same<tmpl::int32_t<3>, tmpl::integral_constant<int32_t, 3>>();
assert_same<tmpl::int64_t<3>, tmpl::integral_constant<int64_t, 3>>();

assert_same<tmpl::uint8_t<3>, tmpl::integral_constant<uint8_t, 3>>();
assert_same<tmpl::uint16_t<3>, tmpl::integral_constant<uint16_t, 3>>();
assert_same<tmpl::uint32_t<3>, tmpl::integral_constant<uint32_t, 3>>();
assert_same<tmpl::uint64_t<3>, tmpl::integral_constant<uint64_t, 3>>();

assert_same<tmpl::size_t<3>, tmpl::integral_constant<size_t, 3>>();
assert_same<tmpl::ptrdiff_t<3>, tmpl::integral_constant<ptrdiff_t, 3>>();
assert_same<tmpl::bool_<true>, tmpl::integral_constant<bool, true>>();
/// [tmpl::integral_constant::abbreviations]

/// [tmpl::list]
static_assert(not std::is_same_v<tmpl::list<Type1, Type2>,
                                 tmpl::list<Type2, Type1>>);
/// [tmpl::list]

/// [tmpl::pair]
assert_same<tmpl::pair<Type1, Type2>::first_type, Type1>();
assert_same<tmpl::pair<Type1, Type2>::second_type, Type2>();
/// [tmpl::pair]

/// [tmpl::type_]
assert_same<tmpl::type_<Type1>::type, Type1>();
/// [tmpl::type_]

// Section: Constants

/// [tmpl::empty_sequence]
assert_same<tmpl::empty_sequence, tmpl::list<>>();
/// [tmpl::empty_sequence]

/// [tmpl::false_type]
assert_same<tmpl::false_type, tmpl::bool_<false>>();
/// [tmpl::false_type]

/// [tmpl::no_such_type_]
assert_same<tmpl::index_of<List1<>, Type1>, tmpl::no_such_type_>();
/// [tmpl::no_such_type_]

/// [tmpl::true_type]
assert_same<tmpl::true_type, tmpl::bool_<true>>();
/// [tmpl::true_type]

// Constructor-like functions for lists

/// [tmpl::filled_list]
assert_same<tmpl::filled_list<Type1, 3, List1>, List1<Type1, Type1, Type1>>();
assert_same<tmpl::filled_list<Type1, 3>, tmpl::list<Type1, Type1, Type1>>();
/// [tmpl::filled_list]

/// [tmpl::integral_list]
assert_same<tmpl::integral_list<int, 3, 2, 1>,
            tmpl::list<tmpl::integral_constant<int, 3>,
                       tmpl::integral_constant<int, 2>,
                       tmpl::integral_constant<int, 1>>>();
/// [tmpl::integral_list]

/// [tmpl::range]
assert_same<tmpl::range<size_t, 4, 7>,
            tmpl::list<tmpl::size_t<4>, tmpl::size_t<5>, tmpl::size_t<6>>>();
assert_same<tmpl::range<size_t, 4, 4>, tmpl::list<>>();
/// [tmpl::range]

/// [tmpl::reverse_range]
assert_same<tmpl::reverse_range<size_t, 7, 4>,
            tmpl::list<tmpl::size_t<7>, tmpl::size_t<6>, tmpl::size_t<5>>>();
assert_same<tmpl::reverse_range<size_t, 7, 7>, tmpl::list<>>();
/// [tmpl::reverse_range]

// Section: Functions for querying lists

/// [tmpl::at_c]
assert_same<tmpl::at_c<List1<Type1, Type2, Type3>, 0>, Type1>();
/// [tmpl::at_c]

/// [tmpl::back]
assert_same<tmpl::back<List1<Type1, Type2, Type3>>, Type3>();
/// [tmpl::back]

/// [tmpl::count_if]
assert_same<tmpl::count_if<List1<Type1, Type2, Type1>,
                           std::is_same<tmpl::_1, Type1>>,
            tmpl::integral_constant<size_t, 2>>();
/// [tmpl::count_if]

/// [tmpl::front]
assert_same<tmpl::front<List1<Type1, Type2, Type3>>, Type1>();
/// [tmpl::front]

/// [tmpl::index_if]
assert_same<tmpl::index_if<List1<Type1, Type2, Type3>,
                           std::is_same<Type3, tmpl::_1>>,
            tmpl::size_t<2>>();
assert_same<tmpl::index_if<List1<Type1, Type3, Type3>,
                           std::is_same<Type3, tmpl::_1>>,
            tmpl::size_t<1>>();
assert_same<tmpl::index_if<List1<Type1>, std::is_same<Type3, tmpl::_1>>,
            tmpl::no_such_type_>();
assert_same<tmpl::index_if<List1<Type1>, std::is_same<Type3, tmpl::_1>, Type2>,
            Type2>();
/// [tmpl::index_if]

/// [tmpl::index_of]
assert_same<tmpl::index_of<List1<Type1, Type2, Type3>, Type3>,
            tmpl::size_t<2>>();
assert_same<tmpl::index_of<List1<Type1, Type3, Type3>, Type3>,
            tmpl::size_t<1>>();
assert_same<tmpl::index_of<List1<Type1>, Type2>,
            tmpl::no_such_type_>();
/// [tmpl::index_of]

/// [tmpl::list_contains]
assert_same<tmpl::list_contains<List1<Type1, Type2>, Type1>,
            tmpl::integral_constant<bool, true>>();
static_assert(tmpl::list_contains_v<List1<Type1, Type2>, Type1>);
static_assert(not tmpl::list_contains_v<List1<Type2, Type2>, Type1>);
/// [tmpl::list_contains]

/// [tmpl::size]
assert_same<tmpl::size<List1<Type1, Type1>>,
            tmpl::integral_constant<unsigned int, 2>>();
/// [tmpl::size]

// Section: Functions producing lists from other lists

{
using lazy_test_arguments =
    tmpl::list<List1<Type1, Type2>, List2<>, List2<Type2>>;
/// [tmpl::append]
assert_same<tmpl::append<List1<Type1, Type2>, List2<>, List2<Type2>>,
            List1<Type1, Type2, Type2>>();
assert_same<tmpl::append<>, tmpl::list<>>();
HAS_LAZY_VERSION(append);
/// [tmpl::append]
}

/// [tmpl::append::bug]
assert_same<tmpl::append<List1<Type1, Type2>, List2<>, Type2>, tmpl::list<>>();
/// [tmpl::append::bug]

/// [tmpl::clear]
assert_same<tmpl::clear<List1<Type1>>, List1<>>();
/// [tmpl::clear]

/// [tmpl::erase_c]
assert_same<tmpl::erase_c<List1<Type1, Type2, Type3>, 1>,
            List1<Type1, Type3>>();
/// [tmpl::erase_c]

{
using lazy_test_arguments = tmpl::list<List1<Type1, Type2, Type1, Type3>,
                                       std::is_same<Type1, tmpl::_1>>;
/// [tmpl::filter]
assert_same<tmpl::filter<List1<Type1, Type2, Type1, Type3>,
                         std::is_same<Type1, tmpl::_1>>,
            List1<Type1, Type1>>();
HAS_LAZY_VERSION(filter);
/// [tmpl::filter]
}

{
using lazy_test_arguments =
    tmpl::list<List1<List1<Type1, List1<Type2>>, List2<List1<Type3>>>>;
/// [tmpl::flatten]
assert_same<
  tmpl::flatten<List1<List1<Type1, List1<Type2>>, List2<List1<Type3>>>>,
  List1<Type1, Type2, List2<List1<Type3>>>>();
HAS_LAZY_VERSION(flatten);
/// [tmpl::flatten]
}

/// [tmpl::join]
assert_same<tmpl::join<List3<List1<Type1, Type2>, List2<>, List2<Type2>>>,
            List1<Type1, Type2, Type2>>();
assert_same<tmpl::join<List1<>>, tmpl::list<>>();
/// [tmpl::join]

/// [tmpl::join::bug]
assert_same<tmpl::join<List3<List1<Type1, Type2>, List2<>, Type2>>,
            tmpl::list<>>();
/// [tmpl::join::bug]

/// [tmpl::join::bug-lazy]
assert_same<
  tmpl::lazy::join<List3<List1<Type1, Type2>, List2<>, List2<Type2>>>::type,
  List3<Type1, Type2, Type2>>();
assert_same<tmpl::lazy::join<List1<>>::type, List1<>>();
/// [tmpl::join::bug-lazy]

/// [tmpl::pop_back]
assert_same<tmpl::pop_back<List1<Type1, Type2, Type3>>, List1<Type1, Type2>>();
assert_same<tmpl::pop_back<List1<Type1, Type2, Type3>,
                           tmpl::integral_constant<unsigned int, 2>>,
            List1<Type1>>();
/// [tmpl::pop_back]

{
using lazy_test_arguments =
    tmpl::list<List1<Type1, Type2, Type3>,
               tmpl::integral_constant<unsigned int, 2>>;
/// [tmpl::pop_front]
assert_same<tmpl::pop_front<List1<Type1, Type2, Type3>>,
            List1<Type2, Type3>>();
assert_same<tmpl::pop_front<List1<Type1, Type2, Type3>,
                            tmpl::integral_constant<unsigned int, 2>>,
            List1<Type3>>();
HAS_LAZY_VERSION(pop_front);
/// [tmpl::pop_front]
}

/// [tmpl::push_back]
assert_same<tmpl::push_back<List1<Type1>, Type2, Type3>,
            List1<Type1, Type2, Type3>>();
/// [tmpl::push_back]

{
using lazy_test_arguments = tmpl::list<List1<Type1>, Type2, Type3>;
/// [tmpl::push_front]
assert_same<tmpl::push_front<List1<Type1>, Type2, Type3>,
            List1<Type2, Type3, Type1>>();
HAS_LAZY_VERSION(push_front);
/// [tmpl::push_front]
}

{
using lazy_test_arguments =
    tmpl::list<List1<Type1, Type2, Type1, Type3>, Type1>;
/// [tmpl::remove]
assert_same<tmpl::remove<List1<Type1, Type2, Type1, Type3>, Type1>,
            List1<Type2, Type3>>();
HAS_LAZY_VERSION(remove);
/// [tmpl::remove]
}

{
using lazy_test_arguments = tmpl::list<List1<Type1, Type2, Type1, Type3>,
                                       std::is_same<Type1, tmpl::_1>>;
/// [tmpl::remove_if]
assert_same<tmpl::remove_if<List1<Type1, Type2, Type1, Type3>,
                            std::is_same<Type1, tmpl::_1>>,
            List1<Type2, Type3>>();
HAS_LAZY_VERSION(remove_if);
/// [tmpl::remove_if]
}

/// [tmpl::replace]
assert_same<tmpl::replace<List1<Type1, Type2, Type1>, Type1, Type3>,
            List1<Type3, Type2, Type3>>();
/// [tmpl::replace]

/// [tmpl::replace:bug]
assert_same<tmpl::apply<tmpl::lazy::replace<tmpl::pin<List1<Type1>>,
                                            tmpl::_1, Type2>,
                        Type1>,
            List1<Type1>>();
assert_same<tmpl::apply<tmpl::bind<tmpl::replace, tmpl::pin<List1<Type1>>,
                                   tmpl::_1, Type2>,
                        Type1>,
            List1<Type2>>();
/// [tmpl::replace:bug]

{
using lazy_test_arguments = tmpl::list<List1<Type1, Type2, Type1>,
                                       std::is_same<Type1, tmpl::_1>, Type3>;
/// [tmpl::replace_if]
assert_same<tmpl::replace_if<List1<Type1, Type2, Type1>,
                             std::is_same<Type1, tmpl::_1>, Type3>,
            List1<Type3, Type2, Type3>>();
HAS_LAZY_VERSION(replace_if);
/// [tmpl::replace_if]
}

{
using lazy_test_arguments = tmpl::list<List1<Type1, Type2, Type3, Type1>>;
/// [tmpl::reverse]
assert_same<tmpl::reverse<List1<Type1, Type2, Type3, Type1>>,
            List1<Type1, Type3, Type2, Type1>>();
HAS_LAZY_VERSION(reverse);
/// [tmpl::reverse]
}

{
using lazy_test_arguments =
    tmpl::list<List1<Type1, Type2, Type3, Type2, Type3, Type3, Type1>, Type3>;
/// [tmpl::split]
assert_same<
  tmpl::split<List1<Type1, Type2, Type3, Type2, Type3, Type3, Type1>, Type3>,
  List1<List1<Type1, Type2>, List1<Type2>, List1<Type1>>>();
HAS_LAZY_VERSION(split);
/// [tmpl::split]
}

{
using lazy_test_arguments =
    tmpl::list<List1<Type1, Type2, Type3>,
               tmpl::integral_constant<unsigned int, 2>>;
/// [tmpl::split_at]
assert_same<tmpl::split_at<List1<Type1, Type2, Type3>,
                           tmpl::integral_constant<unsigned int, 2>>,
            List1<List1<Type1, Type2>, List1<Type3>>>();
HAS_LAZY_VERSION(split_at);
/// [tmpl::split_at]
}

{
using lazy_test_arguments =
    tmpl::list<List2<Type1, Type2, Type3>, List3<Type3, Type2, Type1>,
               make_list1<tmpl::_1, tmpl::_2>>;
/// [tmpl::transform]
assert_same<
  tmpl::transform<List2<Type1, Type2, Type3>, List3<Type3, Type2, Type1>,
                  make_list1<tmpl::_1, tmpl::_2>>,
  List2<List1<Type1, Type3>, List1<Type2, Type2>, List1<Type3, Type1>>>();
HAS_LAZY_VERSION(transform);
/// [tmpl::transform]
}

// Mathematical functions

/// [math_arithmetic]
assert_same<tmpl::plus<tmpl::size_t<10>, tmpl::size_t<3>>::type,
            tmpl::size_t<13>>();
assert_same<tmpl::minus<tmpl::size_t<10>, tmpl::size_t<3>>::type,
            tmpl::size_t<7>>();
assert_same<tmpl::times<tmpl::size_t<10>, tmpl::size_t<3>>::type,
            tmpl::size_t<30>>();
assert_same<tmpl::divides<tmpl::size_t<10>, tmpl::size_t<3>>::type,
            tmpl::size_t<3>>();
assert_same<tmpl::modulo<tmpl::size_t<10>, tmpl::size_t<3>>::type,
            tmpl::size_t<1>>();
assert_same<tmpl::negate<tmpl::int64_t<10>>::type, tmpl::int64_t<-10>>();
/// [math_arithmetic]

/// [math_bitwise]
assert_same<tmpl::complement<tmpl::uint8_t<0b10001111>>::type,
                             tmpl::uint8_t<0b01110000>>();
assert_same<tmpl::bitand_<tmpl::uint8_t<0b00111011>,
                          tmpl::uint8_t<0b01010110>>::type,
                          tmpl::uint8_t<0b00010010>>();
assert_same<tmpl::bitor_<tmpl::uint8_t<0b01100011>,
                         tmpl::uint8_t<0b10100111>>::type,
                         tmpl::uint8_t<0b11100111>>();
assert_same<tmpl::bitxor_<tmpl::uint8_t<0b11000011>,
                          tmpl::uint8_t<0b00000110>>::type,
                          tmpl::uint8_t<0b11000101>>();
assert_same<tmpl::shift_left<tmpl::uint8_t<0b10001110>, tmpl::size_t<3>>::type,
                             tmpl::uint8_t<0b01110000>>();
assert_same<tmpl::shift_right<tmpl::uint8_t<0b10110011>, tmpl::size_t<4>>::type,
                              tmpl::uint8_t<0b00001011>>();
/// [math_bitwise]

/// [math_comparison]
assert_same<tmpl::equal_to<tmpl::size_t<1>, tmpl::size_t<2>>::type,
            tmpl::false_type>();
assert_same<tmpl::not_equal_to<tmpl::size_t<1>, tmpl::size_t<2>>::type,
            tmpl::true_type>();
assert_same<tmpl::greater<tmpl::size_t<1>, tmpl::size_t<2>>::type,
            tmpl::false_type>();
assert_same<tmpl::greater_equal<tmpl::size_t<1>, tmpl::size_t<2>>::type,
            tmpl::false_type>();
assert_same<tmpl::less<tmpl::size_t<1>, tmpl::size_t<2>>::type,
            tmpl::true_type>();
assert_same<tmpl::less_equal<tmpl::size_t<1>, tmpl::size_t<2>>::type,
            tmpl::true_type>();
/// [math_comparison]

/// [math_logical]
assert_same<tmpl::and_<std::true_type, std::false_type>::type,
            tmpl::false_type>();
assert_same<tmpl::or_<std::true_type, std::false_type>::type,
            tmpl::true_type>();
assert_same<tmpl::xor_<std::true_type, std::false_type>::type,
            tmpl::true_type>();
assert_same<tmpl::not_<std::true_type>::type, tmpl::false_type>();
/// [math_logical]

/// [tmpl::identity]
assert_same<tmpl::identity<tmpl::size_t<10>>::type, tmpl::size_t<10>>();
/// [tmpl::identity]

/// [tmpl::max]
assert_same<tmpl::max<tmpl::size_t<10>, tmpl::int32_t<3>>::type,
            tmpl::size_t<10>>();
/// [tmpl::max]

/// [tmpl::min]
assert_same<tmpl::min<tmpl::size_t<10>, tmpl::size_t<3>>::type,
            tmpl::size_t<3>>();
/// [tmpl::min]

/// [tmpl::next]
assert_same<tmpl::next<tmpl::size_t<10>>::type, tmpl::size_t<11>>();
/// [tmpl::next]

/// [tmpl::prev]
assert_same<tmpl::prev<tmpl::size_t<10>>::type, tmpl::size_t<9>>();
/// [tmpl::prev]

// Miscellaneous functions

/// [tmpl::always]
assert_same<tmpl::always<Type1>::type, Type1>();
/// [tmpl::always]

/// [tmpl::apply]
assert_same<tmpl::apply<std::is_convertible<tmpl::_1, tmpl::_2>,
                        const char*, std::string>,
            std::true_type>();
assert_same<tmpl::apply<std::is_convertible<tmpl::_2, tmpl::_1>,
                        const char*, std::string>,
            std::false_type>();
/// [tmpl::apply]

/// [tmpl::count]
assert_same<tmpl::count<Type1, Type2, Type1>,
            tmpl::integral_constant<unsigned int, 3>>();
/// [tmpl::count]

/// [tmpl::if_]
assert_same<tmpl::if_<std::true_type, Type1, Type2>::type, Type1>();
/// [tmpl::if_]

/// [tmpl::if_c]
assert_same<tmpl::if_c<true, Type1, Type2>::type, Type1>();
/// [tmpl::if_c]

/// [tmpl::is_set]
assert_same<tmpl::is_set<Type1, Type2, Type3>, tmpl::true_type>();
assert_same<tmpl::is_set<Type1, Type2, Type1>, tmpl::false_type>();
assert_same<tmpl::is_set<>, tmpl::true_type>();
/// [tmpl::is_set]

/// [tmpl::repeat]
assert_same<tmpl::repeat<Wrapper, tmpl::size_t<3>, Type1>,
            Wrapper<Wrapper<Wrapper<Type1>>>>();
assert_same<tmpl::repeat<Wrapper, tmpl::size_t<0>, Type1>, Type1>();
/// [tmpl::repeat]

/// [tmpl::repeat:lazy]
assert_same<tmpl::lazy::repeat<Wrapper, tmpl::size_t<3>, Type1>::type,
            Wrapper<Wrapper<Wrapper<Type1>>>>();
/// [tmpl::repeat:lazy]

/// [tmpl::sizeof_]
assert_same<tmpl::sizeof_<double>::type,
            tmpl::integral_constant<unsigned int, sizeof(double)>>();
/// [tmpl::sizeof_]

/// [tmpl::substitute]
assert_same<tmpl::substitute<List1<List2<tmpl::_1, tmpl::_2, tmpl::_3>,
                                   tmpl::args<0>, tmpl::args<1>, tmpl::args<2>>,
                             List3<Type1, Type2, Type3>>,
            List1<List2<tmpl::_1, tmpl::_2, Type3>, Type1, Type2, Type3>>();
/// [tmpl::substitute]

/// [tmpl::type_from]
assert_same<tmpl::type_from<tmpl::type_<Type1>>, Type1>();
/// [tmpl::type_from]

/// [tmpl::wrap]
assert_same<tmpl::wrap<List1<Type1, Type2>, List2>, List2<Type1, Type2>>();
/// [tmpl::wrap]

/// [tmpl::wrap:lazy]
assert_same<tmpl::lazy::wrap<List1<Type1, Type2>, List2>::type,
            List2<Type1, Type2>>();
/// [tmpl::wrap:lazy]
}
} // namespace
// clang-format on
