\cond NEVER
Distributed under the MIT License.
See LICENSE.txt for details.
\endcond
# Metaprogramming with Brigand {#brigand}

\tableofcontents

[comment]: # (The \pars improve the spacing in the generated document when)
[comment]: # (many \snippets are involved.)

\par
In SpECTRE, most complex TMP is done using the [Brigand metaprogramming
library](https://github.com/edouarda/brigand).  Brigand is a collection of
templated classes, type aliases, and functions, primarily intended to help with
the manipulation and use of lists of types.  This document is organized to
roughly parallel the structure of the C++ standard, rather than following
Brigand's classifications.

\par
Brigand provides all of its functionality in the `brigand` namespace, but in
SpECTRE we have aliased this namespace to `tmpl`, and the latter should be
preferred.

\par
All functionality described here is provided by SpECTRE's Brigand wrapper
header:
\snippet Test_TMPLDocumentation.cpp include

\par
Examples in this document use the following declarations and definitions:
\snippet Test_TMPLDocumentation.cpp example_declarations


\section Containers

\par
Brigand provides container classes with the sole purpose of wrapping other
things.


\subsection integral_constant

\par
Very similar to std::integral_constant, except that "for maximum C++11
compatibility" all the `constexpr` specifiers have been omitted.
\snippet Test_TMPLDocumentation.cpp tmpl::integral_constant

\par
Brigand supplies type aliases for constants of some specific types:
\snippet Test_TMPLDocumentation.cpp tmpl::integral_constant::abbreviations


\subsection integral_list

\par
Shorthand for a \ref list of \ref integral_constant ""s of the same type
\snippet Test_TMPLDocumentation.cpp tmpl::integral_list


\subsection list

\par
An empty struct templated on a parameter pack, with no additional
functionality.
\snippet Test_TMPLDocumentation.cpp tmpl::list

\par
Most metafunctions that operate on lists will work on any struct template.


\subsection type_

\par
A struct templated on a single type `T` containing an alias `type` to `T`.
\snippet Test_TMPLDocumentation.cpp tmpl::type_

\par
When extracting the type, programmers are encouraged to use \ref type_from.


\section Constants

\par
Brigand defines a few concrete types and type aliases.


\subsection empty_sequence

\par
An empty \ref list.
\snippet Test_TMPLDocumentation.cpp tmpl::empty_sequence

\par
Prefer just writing an empty list.


\subsection false_type

\par
An \ref integral_constant representing `false`.  Similar to std::false_type.
\snippet Test_TMPLDocumentation.cpp tmpl::false_type


\subsection no_such_type_

\par
An empty struct returned as the failure case for various searching operations.
\snippet Test_TMPLDocumentation.cpp tmpl::no_such_type_


\subsection true_type

\par
An \ref integral_constant representing `true`.  Similar to std::true_type.
\snippet Test_TMPLDocumentation.cpp tmpl::true_type


\section list_query Functions for querying lists

\par
These tend to be similar to const member functions in the STL and the
non-modifying sequence operations in `<algorithm>`.  They are most frequently
used with \ref list, but similar classes will also work.


\subsection at_c

\par
Retrieves a given element of a list, similar to operator[] of the STL
containers.  The index number is supplied as an `unsigned int`.
\snippet Test_TMPLDocumentation.cpp tmpl::at_c


\subsection back

\par
Retrieves the last element of a list.
\snippet Test_TMPLDocumentation.cpp tmpl::back


\subsection front

\par
Retrieves the first element of a list.
\snippet Test_TMPLDocumentation.cpp tmpl::front


\subsection index_of

\par
Finds the index as a size_t \ref integral_constant of the first occurrence of a
type in a list.  Returns \ref no_such_type_ if the type is not found.
\snippet Test_TMPLDocumentation.cpp tmpl::index_of


\subsection list_contains

\par
Checks whether a particular type is contained in a list, returning an
integral_constant of bool.
\snippet Test_TMPLDocumentation.cpp tmpl::list_contains

\note
This is not a Brigand metafunction.  It is implemented in SpECTRE.


\subsection size

\par
Returns the size of a list as an \ref integral_constant of type `unsigned int`.
\snippet Test_TMPLDocumentation.cpp tmpl::size


\section list_to_list Functions producing lists from other lists

\par
These tend to be similar to non-const member functions in the STL and the
mutating sequence operations in `<algorithm>`, but due to the nature of
metaprogramming all return a new list rather than modifying an argument.  They
are most frequently used with \ref list, but similar classes will also work.


\subsection append

\par
Appends the contents of several lists to the contents of a list.
\snippet Test_TMPLDocumentation.cpp tmpl::append

\par
For a version taking its arguments as a list of lists, see \ref join.

\warning
A flaw in the implementation makes use of this metafunction very error-prone.
Prefer \ref push_back or \ref push_front when possible.
\snippet Test_TMPLDocumentation.cpp tmpl::append::bug


\subsection clear

\par
Produces list with the same head but no elements.
\snippet Test_TMPLDocumentation.cpp tmpl::clear


\subsection erase_c

\par
Produces a copy of a list with the element at the given index (passed as an
`unsigned int`) removed.
\snippet Test_TMPLDocumentation.cpp tmpl::erase_c


\subsection join

\par
Appends to a list in the same manner as \ref append, but takes a list of lists
instead of multiple arguments.
\snippet Test_TMPLDocumentation.cpp tmpl::join

\warning
A flaw in the implementation makes use of this metafunction very error-prone.
Prefer \ref push_back or \ref push_front when possible.
\snippet Test_TMPLDocumentation.cpp tmpl::join::bug


\subsection pop_back

\par
Remove elements from the end of a list.  The number of elements to remove is
supplied as an \ref integral_constant and defaults to 1.
\snippet Test_TMPLDocumentation.cpp tmpl::pop_back


\subsection pop_front

\par
Remove elements from the beginning of a list.  The number of elements to remove
is supplied as an \ref integral_constant and defaults to 1.
\snippet Test_TMPLDocumentation.cpp tmpl::pop_front


\subsection push_back

\par
Appends types to a list.
\snippet Test_TMPLDocumentation.cpp tmpl::push_back


\subsection push_front

\par
Prepends types to a list.  The order of the prepended items is retained: they
are pushed as a unit, not one-by-one.
\snippet Test_TMPLDocumentation.cpp tmpl::push_front


\subsection split_at

\par
Given a list and an integer \f$N\f$ (supplied as an \ref integral_constant),
returns a list of two of lists, the first containing the first \f$N\f$
elements, and the second containing the remaining elements.
\snippet Test_TMPLDocumentation.cpp tmpl::split_at


\section TODO

\subsection type_from


```
    182 tmpl::_1
    157 tmpl::for_each
    155 tmpl::conditional_t
    103 tmpl::transform
     98 tmpl::bind
     94 tmpl::flatten
     77 tmpl::type_from
     63 tmpl::pin
     42 tmpl::remove_duplicates
     30 tmpl::filter
     17 tmpl::list_difference
     14 tmpl::at
      8 tmpl::position_of_first_v
      8 tmpl::next
      8 tmpl::all
      8 tmpl::_state
      8 tmpl::_element
      7 tmpl::fold
      6 tmpl::make_sequence
      6 tmpl::erase
      6 tmpl::count_if
      5 tmpl::sort
      5 tmpl::parent
      5 tmpl::max
      5 tmpl::found
      4 tmpl::remove_if
      4 tmpl::equal_members
      3 tmpl::wrap
      3 tmpl::replace_at
      3 tmpl::position_of_first
      3 tmpl::not_
      3 tmpl::map
      3 tmpl::has_key
      3 tmpl::and_
      2 tmpl::plus
      2 tmpl::pair
      2 tmpl::is_set
      2 tmpl::get_source
      2 tmpl::get_destination
      2 tmpl::edge
      2 tmpl::apply
      2 tmpl::any
      2 tmpl::_2
      1 tmpl::reverse
      1 tmpl::replace
      1 tmpl::remove
      1 tmpl::range
      1 tmpl::or_
      1 tmpl::less
      1 tmpl::lazy
      1 tmpl::insert
      1 tmpl::index_if
      1 tmpl::find
      1 tmpl::filled_list
      1 tmpl::enumerated_fold
      1 tmpl::count
      1 tmpl::_3

as_tuple
```
