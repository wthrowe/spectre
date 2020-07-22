\cond NEVER
Distributed under the MIT License.
See LICENSE.txt for details.
\endcond
# Metaprogramming with Brigand {#brigand}

\note
This document covers Brigand as of commit
[85baf9e685eb0c942764b7224fa1ce034bb3beba](https://github.com/edouarda/brigand/commit/85baf9e685eb0c942764b7224fa1ce034bb3beba)
in Summer 2017.  There have been only minor changes since then.

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


\section Metafunctions

\par
A metafunction is an analog of a familiar C++ function that is coded in the C++
type system.  It turns out that, using metafunction programming, it is possible
to perform arbitrary computations at compile time.

\par
There are multiple ways to encode a calculation in the type system.  When using
Brigand, the relevant ones are eager and lazy metafunctions.


\subsection lazy Eager and lazy metafunctions

\par
Metafunctions commonly appear in two forms: eager and lazy.  Lazy functions are
templated structs (or templated aliases to structs) with a `type` member
function that indicates the result:
\snippet Test_TMPLDocumentation.cpp metafunctions:lazy
The type traits in the standard library, such as std::is_same, are lazy
metafunctions.

\par
Eager metafunctions are aliases to their result types.  As a trivial case,
struct templates can be viewed as eager metafunctions returning themselves.  An
eager version of the previous example could be implemented as:
\snippet Test_TMPLDocumentation.cpp metafunctions:eager
The standard library provides eager versions of some of its metafunctions
(generally those that modify a type, rather than predicates) using an `_t`
suffix.  When both versions are provided, it is often convenient (and less
error prone!) to define the eager version in terms of the lazy version:
\snippet Test_TMPLDocumentation.cpp metafunctions:eager_from_lazy
And the two definitions agree:
\snippet Test_TMPLDocumentation.cpp metafunctions:agreement

\note
The standard library also provides versions of many of its type traits with an
`_v` suffix.  These evaluate to compile-time *values*, rather than types.  They
can be useful for metaprogramming, but are not the types of metafunctions being
discussed here.

\par
Eager metafunctions are usually more convenient to use, so what is the point of
additionally creating lazy ones?  The answer is that lazy metafunctions can be
used as compile-time functors.  As a simple example, we can write a function
that calls an arbitrary lazy metafunction several times
\snippet Test_TMPLDocumentation.cpp metafunctions:lazy_call_metafunction
and get the expected output:
\snippet Test_TMPLDocumentation.cpp metafunctions:lazy_call_metafunction_assert
A similar attempt with an eager metafunction fails, because the function is
evaluated too early, acting as a function composition, rather than a lambda:
\snippet Test_TMPLDocumentation.cpp metafunctions:eager_call_metafunction
\snippet Test_TMPLDocumentation.cpp metafunctions:eager_call_metafunction_assert
(In this simple case we could have used a template template parameter to pass
the eager metafunction in a form more similar to a runtime lambda, but the
possibilities for generic manipulation of parameter lists containing template
template parameters are limited, so their use must be minimized in complex
metaprogramming.)

\note
In practice, lazy metafunctions are often implemented as empty structs
inheriting from other lazy metafunctions.  The entire inheritance chain then
inherits a `type` alias from the ultimate base class.

\par
Most of the standard Brigand functions are eager, but many have lazy versions
in the nested `lazy` namespace.  Thesw are indicated by calls to the
`HAS_LAZY_VERSION` macro in the examples below.


\subsection metalambdas Brigand metalambdas

\par
This use of lazy metafunctions is too limited for general use, however, because
it requires the definition of a new templated struct for every new function.
Brigand uses a more general notation, known as metalambdas.  A metalambda is a
(possibly nested set of) lazy metafunctions with some template arguments
replaced by the placeholders `_1`, `_2`, etc.  These are the first, second,
etc., arguments of the metalambda, and will be replaced by the actual arguments
when the lambda is used.  The lazy nature of the metafunctions prevents them
from prematurely evaluating to results based on the literal placeholder types.
The \ref apply function can be used to evaluate a metalambda with specified
arguments, any many other Brigand functions take metalambdas that are evaluated
internally.


\subsection metalambda_structure Evaluation of metalambdas

\note
None of the terminology introduced in this section is standard.

\par
A metalambda is always evaluated in some *context* describing the possible
argument substitutions to be performed.  The context is a list of arguments
and, possibly, a parent context.  When evaluating a metalambda, for example by
calling \ref apply, the initial context contains the passed set of arguments
with no parent context.

\par
There are eight forms that a metalambda can take: an argument, a lazy
expression, a bind expression, a pin expression, a defer expression, a parent
expression, a constant, or a metaclosure.

\subsubsection args Argument

\par
An argument is one of the structs `_1`, `_2`, or `args<N>` for `unsigned int`
N.  The additional aliases `_3`, `_4`, ..., `_9` are provided to `args<2>`,
`args<3>`, ..., `args<8>`.
\snippet Test_TMPLDocumentation.cpp tmpl::args
When evaluated, they give the first (`_1`), second (`_2`), or zero-indexed Nth
(`args<N>`) element of the context's argument list.
\snippet Test_TMPLDocumentation.cpp tmpl::args:eval
Additionally, `_state` and `_element` are aliased to `_1` and `_2`, primarily
for use with \ref fold.

\par
If there are too few arguments in the argument list, `_1` and `_2` remain
unevaluated, but `args<N>` (including `args<0>` and `args<1>`) gives an error.
Do not rely on this behavior.

\subsubsection metalambda_lazy Lazy expression

\par
A lazy expression is a fully-specialized struct template with only type
template parameters that is not a specialization of \ref pin "pin", \ref defer
"defer", or \ref parent "parent" and is not a \ref metalambda_metaclosure
"metaclosure".  When evaluated, each of its template parameters is evaluated as
a metafunction and replaced by the result, and then the struct's `type` type
alias is the result of the full lazy-expression.
\snippet Test_TMPLDocumentation.cpp metalambda_lazy

\subsubsection bind Bind expression

\par
A bind expression is a specialization of `bind`.  It wraps an eager
metafunction and its arguments.  When evaluated, the arguments are each
evaluated as metafunctions, and then the results are passed to the eager
metafunction.
\snippet Test_TMPLDocumentation.cpp tmpl::bind

\note
The `bind` metafunction does not convert an eager metafunction to a lazy one.
It is handled specially in the evaluation code.

\subsubsection pin Pin expression

\par
A pin expression is a specialization of `pin`.  Evaluating a pin expression
gives the argument to `pin`.  This can be used to force a type to be treated as
a \ref metalambda_constant "constant", even if it would normally be treated as
a different type of metalambda (usually a \ref metalambda_lazy
"lazy expression").
\snippet Test_TMPLDocumentation.cpp tmpl::pin

\subsubsection defer Defer expression

\par
A defer expression is a specialization of `defer`.  It does not evaluate its
argument, but results in a \ref metalambda_metaclosure "metaclosure" containing
the passed metalambda and the current evaluation context.
\snippet Test_TMPLDocumentation.cpp tmpl::defer

\par
The primary purposes for `defer` are constructing metalambdas to pass to other
metafunctions and preventing "speculative" evaluation of a portion of a
metalambda that is not valid for some arguments.  See the examples below, in
particular \ref multiplication_table, \ref maybe_first, and \ref
column_with_zeros.

\warning
The metalambda contained in a `defer` must be a \ref metalambda_lazy
"lazy expression" or a \ref bind "bind expression".  This is presumably a bug.
If another type is needed, it can be wrapped in \ref always.

\subsubsection parent Parent expression

\par
A parent expression is a specialization of `parent`.  It evaluates its
argument, replacing the current context with its parent.  This provides access
to the captured arguments in a metaclosure.
\snippet Test_TMPLDocumentation.cpp tmpl::parent

\warning
Do not call `parent` outside of a metaclosure context.  This results in an
empty evaluation context, causing unintuitive changes to the evaluation rules.
(Most, but not all, expressions are left unevaluated in such a context.)  Use
\ref pin "pin" to prevent evaluation.

\par

\warning
There is a bug that prevents `parent` from working in a metaclosure being
evaluated in a metaclosure context.  In some cases this can be worked around by
evaluating the metaclosure in the parent of the metaclosure context.
FIXME add snippet for bug

\subsubsection metalambda_constant Constant

\par
A constant metalambda is any type that is not a struct template with only type
template parameters, a specialization of \ref bind "bind", or a metaclosure.  A
constant metalambda evaluates to itself.
\snippet Test_TMPLDocumentation.cpp metalambda_constant

\subsubsection metalambda_metaclosure Metaclosure

\par
A metaclosure is an opaque type produced by \ref defer "defer", containing a
metalambda and an evaluation context.  When a metaclosure is evaluated, it
evaluates the packaged metalambda in the current evaluation context with the
parent context replaced by the packaged context.  See \ref defer and \ref
parent for examples.


\subsection Examples


\subsubsection evens

\par
Finds all numbers in a list that are even.
\snippet Test_TMPLDocumentation.cpp metafunctions:evens
\snippet Test_TMPLDocumentation.cpp metafunctions:evens:asserts


\subsubsection maybe_first

\par
Returns the first element of a list, or \ref no_such_type_ if the list is
empty.
\snippet Test_TMPLDocumentation.cpp metafunctions:maybe_first
\snippet Test_TMPLDocumentation.cpp metafunctions:maybe_first:asserts

\par
This example demonstrates the use of \ref defer "defer" to lazily evaluate a
branch of the \ref if_, preventing an attempt to evaluate
`tmpl::front<tmpl::list<>>`.


\subsubsection factorial

\par
Calculates the factorial using a simple metalambda passed to a \ref fold.
\snippet Test_TMPLDocumentation.cpp metafunctions:factorial
\snippet Test_TMPLDocumentation.cpp metafunctions:factorial:asserts


\subsubsection multiplication_table

\par
Constructs a multiplication table.
\snippet Test_TMPLDocumentation.cpp metafunctions:multiplication_table
\snippet Test_TMPLDocumentation.cpp metafunctions:multiplication_table:asserts

\par
This demonstrates the use of \ref defer "defer" to pass a closure as an
argument to a metafunction (\ref transform "tmpl::lazy::transform"), while
capturing an argument from the outer context (the metalambda evaluated for the
outer `tmpl::transform`).


\subsubsection column_with_zeros

\par
Extracts a column from a row-major matrix, extending any short rows with zeros.
\snippet Test_TMPLDocumentation.cpp metafunctions:column_with_zeros
\snippet Test_TMPLDocumentation.cpp metafunctions:column_with_zeros:asserts

\par
This example shows another use of \ref defer "defer" to avoid evaluating an
invalid expression, similar to \ref maybe_first.  The use of an \ref args
"argument" in the deferred branch makes this case more complicated: a \ref
parent "parent" expression is used to access the context where the \ref defer
"defer" occurs to avoid having to pass the argument explicitly using the \ref
apply call.

\par
This is the "apply-defer-parent" pattern for lazy evaluation.  A \ref parent
"parent" is placed immediately inside a \ref defer "defer" (with an \ref always
to work around a Brigand bug) with a (not immediately) surrounding \ref apply.
The \ref parent "parent" causes its contents to be executed (when called by
\ref apply "apply") in the context where the \ref defer "defer" was evaluated,
so the deferral is unobservable by the contents.


\subsubsection factorial_recursion

\par
Again calculates the factorial, but using a recursive algorithm.  After some
setup code to start the recursion, the recursive metalambda is called with
itself as the first argument (as a plain lambda, not a closure).  The other two
arguments are an accumulator and the number of remaining iterations.
\snippet Test_TMPLDocumentation.cpp metafunctions:factorial_recursion
\snippet Test_TMPLDocumentation.cpp metafunctions:factorial_recursion:asserts

\par
This again uses the "apply-defer-parent" pattern to prevent "speculative"
evaluation of conditional branches.  In this example, speculative evaluation of
the branch is invalid because it would recurse infinitely.


\subsubsection primes

\par
Generates a list of prime numbers using the sieve of Eratosthenes.  This
example defines three helper metafunctions.  Two, `zero` and `replace_at`, are
defined only for clarity's sake and could be inlined.  The third,
`range_from_types`, is not easily inlinable, and works around Brigand's lack of
sequence generating functions without non-type template parameters.
\snippet Test_TMPLDocumentation.cpp metafunctions:primes
\snippet Test_TMPLDocumentation.cpp metafunctions:primes:asserts


\section function_docs Brigand functions

\subsection Containers

\par
Brigand provides container classes with the sole purpose of wrapping other
things.


\subsubsection integral_constant

\par
Very similar to std::integral_constant, except that the `constexpr` specifiers
on the member functions have been omitted.
\snippet Test_TMPLDocumentation.cpp tmpl::integral_constant

\par
Brigand supplies type aliases for constants of some specific types:
\snippet Test_TMPLDocumentation.cpp tmpl::integral_constant::abbreviations

\par
Most metafunctions that accept integral_constants will accept any type with a
`value` static member variable.

\par
Because of the `type` type alias, integral_constants behave like lazy
metafunctions returning themselves.  Most lazy metafunctions producing an
integral_constant will actually inherit from their result, so `value` will be
directly available without needing to go through the `type` alias.

\remark
Prefer std::integral_constant, except for the convenience wrapper
`tmpl::size_t` or when necessary for type equality comparison.


\subsubsection list

\par
An empty struct templated on a parameter pack, with no additional
functionality.
\snippet Test_TMPLDocumentation.cpp tmpl::list

\par
Most metafunctions that operate on lists will work on any struct template.


\subsubsection pair

\par
A pair of types, with easy access to each type in the pair.
\snippet Test_TMPLDocumentation.cpp tmpl::pair


\subsubsection type_

\par
A struct templated on a single type `T` containing an alias `type` to `T`.
\snippet Test_TMPLDocumentation.cpp tmpl::type_

\par
When extracting the type, programmers are encouraged to use \ref type_from.


\subsection Constants

\par
Brigand defines a few concrete types and type aliases.


\subsubsection empty_sequence

\par
An empty \ref list.
\snippet Test_TMPLDocumentation.cpp tmpl::empty_sequence

\remark
Prefer just writing an empty list.


\subsubsection false_type

\par
An \ref integral_constant representing `false`.  Similar to std::false_type.
\snippet Test_TMPLDocumentation.cpp tmpl::false_type

\remark
Prefer std::false_type.


\subsubsection no_such_type_

\par
An empty struct returned as the failure case for various searching operations.
\snippet Test_TMPLDocumentation.cpp tmpl::no_such_type_


\subsubsection true_type

\par
An \ref integral_constant representing `true`.  Similar to std::true_type.
\snippet Test_TMPLDocumentation.cpp tmpl::true_type

\remark
Prefer std::true_type.


\subsection list_constructor Constructor-like functions for lists

\par
These functions produce \ref list ""s from non-list values.  They are often
similar to constructors in the STL.


\subsubsection filled_list

\par
Creates a list containing a given number (passed as an `unsigned int`) of the
same type.  The head of the list defaults to \ref list.
\snippet Test_TMPLDocumentation.cpp tmpl::filled_list


\subsubsection integral_list

\par
Shorthand for a \ref list of \ref integral_constant ""s of the same type
\snippet Test_TMPLDocumentation.cpp tmpl::integral_list

\remark
Prefer std::integer_sequence when used for pack expansion.  Prefer this when
the contents need to be manipulated for more complicated metaprogramming.


\subsubsection range

\par
Produces a \ref list of \ref integral_constant ""s representing adjacent
ascending integers, including the starting value and excluding the ending
value.
\snippet Test_TMPLDocumentation.cpp tmpl::range

\see \ref reverse_range


\subsubsection reverse_range

\par
Produces a \ref list of \ref integral_constant ""s representing adjacent
descending integers, including the starting value and excluding the ending
value.
\snippet Test_TMPLDocumentation.cpp tmpl::reverse_range

\see \ref range


\subsection list_query Functions for querying lists

\par
These tend to be similar to const member functions in the STL and the
non-modifying sequence operations in `<algorithm>`.  They are most frequently
used with \ref list, but similar classes will also work.


\subsubsection at_c

\par
Retrieves a given element of a list, similar to operator[] of the STL
containers.  The index number is supplied as an `unsigned int`.
\snippet Test_TMPLDocumentation.cpp tmpl::at_c


\subsubsection back

\par
Retrieves the last element of a list.
\snippet Test_TMPLDocumentation.cpp tmpl::back


\subsubsection count_if

\par
Returns the number of elements of a list satisfying a predicate.
\snippet Test_TMPLDocumentation.cpp tmpl::count_if

\note
If the predicate is neither a \ref bind "bind expression" with a single
argument of `_1` nor a \ref metalambda_lazy "lazy expression" with a single
argument of `_1`, then a bug causes the `type` of the result of the predicate
to be used instead of the result itself.  As long as the predicate returns an
\ref integral_constant or a std::integral_constant this does not matter, as
`type` is a no-op for those classes.
\snippet Test_TMPLDocumentation.cpp tmpl::count_if:bug:definitions
\snippet Test_TMPLDocumentation.cpp tmpl::count_if:bug:asserts


\subsubsection front

\par
Retrieves the first element of a list.
\snippet Test_TMPLDocumentation.cpp tmpl::front


\subsubsection index_if

\par
Finds the index as a `size_t` \ref integral_constant of the first type in a
list satisfying a \ref metalambdas "predicate".  Returns a supplied type,
defaulting to \ref no_such_type_ if no elements match.
\snippet Test_TMPLDocumentation.cpp tmpl::index_if


\subsubsection index_of

\par
Finds the index as a `size_t` \ref integral_constant of the first occurrence of
a type in a list.  Returns \ref no_such_type_ if the type is not found.
\snippet Test_TMPLDocumentation.cpp tmpl::index_of


\subsubsection list_contains

\par
Checks whether a particular type is contained in a list, returning an
\ref integral_constant of `bool`.
\snippet Test_TMPLDocumentation.cpp tmpl::list_contains

\note
This is not a Brigand metafunction.  It is implemented in SpECTRE.


\subsubsection size

\par
Returns the size of a list as an \ref integral_constant of type `unsigned int`.
\snippet Test_TMPLDocumentation.cpp tmpl::size


\subsection list_to_list Functions producing lists from other lists

\par
These tend to be similar to non-const member functions in the STL and the
mutating sequence operations in `<algorithm>`, but due to the nature of
metaprogramming all return a new list rather than modifying an argument.  They
are most frequently used with \ref list, but similar classes will also work.


\subsubsection append

\par
Appends the contents of several lists to the contents of a list.
\snippet Test_TMPLDocumentation.cpp tmpl::append

\par
For a version taking its arguments as a list of lists, see \ref join.

\warning
A flaw in the implementation makes use of this metafunction very error-prone.
Prefer \ref push_back or \ref push_front when possible.
\snippet Test_TMPLDocumentation.cpp tmpl::append::bug


\subsubsection clear

\par
Produces list with the same head but no elements.
\snippet Test_TMPLDocumentation.cpp tmpl::clear

\remark
If the head is known, prefer writing it explicitly.  If the head is irrelevant,
write an empty \ref list.


\subsubsection erase_c

\par
Produces a copy of a list with the element at the given index (passed as an
`unsigned int`) removed.
\snippet Test_TMPLDocumentation.cpp tmpl::erase_c


\subsubsection filter

\par
Removes all types not matching a \ref metalambdas "predicate" from a list.
\snippet Test_TMPLDocumentation.cpp tmpl::filter

\see \ref remove_if


\subsubsection flatten

\par
Given a list, recursively inlines the contents of elements that are lists with
the same head.
\snippet Test_TMPLDocumentation.cpp tmpl::flatten


\subsubsection join

\par
Appends to a list in the same manner as \ref append, but takes a list of lists
instead of multiple arguments.
\snippet Test_TMPLDocumentation.cpp tmpl::join

\warning
A flaw in the implementation makes use of this metafunction very error-prone.
Prefer \ref push_back or \ref push_front when possible.
\snippet Test_TMPLDocumentation.cpp tmpl::join::bug

\par

\warning
This metafunction has a lazy version, but its behavior does not match the eager
version, as it determines the head of the resulting list differently.
\snippet Test_TMPLDocumentation.cpp tmpl::join::bug-lazy


\subsubsection pop_back

\par
Remove elements from the end of a list.  The number of elements to remove is
supplied as an \ref integral_constant and defaults to 1.
\snippet Test_TMPLDocumentation.cpp tmpl::pop_back


\subsubsection pop_front

\par
Remove elements from the beginning of a list.  The number of elements to remove
is supplied as an \ref integral_constant and defaults to 1.
\snippet Test_TMPLDocumentation.cpp tmpl::pop_front


\subsubsection push_back

\par
Appends types to a list.
\snippet Test_TMPLDocumentation.cpp tmpl::push_back


\subsubsection push_front

\par
Prepends types to a list.  The order of the prepended items is retained: they
are pushed as a unit, not one-by-one.
\snippet Test_TMPLDocumentation.cpp tmpl::push_front


\subsubsection remove

\par
Removes all occurrences of a given type from a list.
\snippet Test_TMPLDocumentation.cpp tmpl::remove


\subsubsection remove_if

\par
Removes all types matching a \ref metalambdas "predicate" from a list.
\snippet Test_TMPLDocumentation.cpp tmpl::remove_if

\see \ref filter


\subsubsection replace

\par
Replaces all occurrences of one type in a list with another.
\snippet Test_TMPLDocumentation.cpp tmpl::replace

\warning
This metafunction has a lazy version, but it is broken because it is
implemented as a non-trivial type alias.  Use \ref bind "bind" on the eager
version instead.
\snippet Test_TMPLDocumentation.cpp tmpl::replace:bug


\subsubsection replace_if

\par
Replaces all types matching a \ref metalambdas "predicate" in a list with a
given type.
\snippet Test_TMPLDocumentation.cpp tmpl::replace_if


\subsubsection reverse

\par
Reverses the order of types in a list.
\snippet Test_TMPLDocumentation.cpp tmpl::reverse


\subsubsection split

\par
Splits a list into parts separated by a specified delimiter, discarding empty
parts.
\snippet Test_TMPLDocumentation.cpp tmpl::split


\subsubsection split_at

\par
Given a list and an integer \f$N\f$ (supplied as an \ref integral_constant),
returns a list of two of lists, the first containing the first \f$N\f$
elements, and the second containing the remaining elements.
\snippet Test_TMPLDocumentation.cpp tmpl::split_at


\subsubsection transform

\par
Given a list, calls a \ref metalambdas "metalambda" on each type in the list,
collecting the results in a new list.  If additional lists are supplied,
elements from those lists are passed as additional arguments to the metalambda.
\snippet Test_TMPLDocumentation.cpp tmpl::transform


\subsection math Mathematical functions

\par
These perform the same operations at their language counterparts, but on \ref
integral_constant ""s (or anything else with a `value` static member type of
type `value_type`).  The results inherit from \ref integral_constant ""s of
types noted below.

\par
These are all lazy metafunctions.


\subsubsection math_arithmetic Arithmetic operators

\par
These operations return classes inheriting from \ref integral_constant ""s of
the same type as the result of the language operator on their arguments.  The
integral promotion and conversion rules are applied.  (Contrast the \ref
math_bitwise "bitwise operators".)
\snippet Test_TMPLDocumentation.cpp math_arithmetic

\par
The standard library runtime functors have the same names for std::plus,
std::minus, std::divides, and std::negate, but the other two are
std::multiplies and std::modulus.


\subsubsection math_bitwise Bitwise operators

\par
These operations return classes inheriting from \ref integral_constant ""s of
the same type as their first argument's `value`.  This is *not* generally the
same type as the language operator, even when the types of the values of both
arguments are the same.  (The integer promotion and conversion rules are not
applied.)
\snippet Test_TMPLDocumentation.cpp math_bitwise

\par
The standard library runtime functors are called std::bit_not, std::bit_and,
std::bit_or, and std::bit_xor.


\subsubsection math_comparison Comparison operators

\par
These operations return classes inheriting from \ref integral_constant ""s of
`bool`.
\snippet Test_TMPLDocumentation.cpp math_comparison

\par
The standard library runtime functors have the same names, such as
std::equal_to.


\subsubsection math_logical Logical operators

\par
These operations return classes inheriting from \ref integral_constant ""s of
`bool`.  They should only be used on types wrapping `bool`s.
\snippet Test_TMPLDocumentation.cpp math_logical

\par
The standard library runtime functors are called std::logical_and,
std::logical_or, and std::logical_not.  The xor operation is equivalent to \ref
math_comparison "not_equal_to".


\subsubsection identity

\par
The identity function.  Unlike most math functions, this returns the same type
as its argument, even if that is not an \ref integral_constant.
\snippet Test_TMPLDocumentation.cpp tmpl::identity

\see \ref always


\subsubsection max

\par
Computes the larger of two values, returning an \ref integral_constant of the
common type of its arguments.
\snippet Test_TMPLDocumentation.cpp tmpl::max


\subsubsection min

\par
Computes the smaller of two values, returning an \ref integral_constant of the
common type of its arguments.
\snippet Test_TMPLDocumentation.cpp tmpl::min


\subsubsection next

\par
Computes the passed value plus one, returning an \ref integral_constant of the
same type as its argument.
\snippet Test_TMPLDocumentation.cpp tmpl::next


\subsubsection prev

\par
Computes the passed value minus one, returning an \ref integral_constant of the
same type as its argument.
\snippet Test_TMPLDocumentation.cpp tmpl::prev


\subsection misc Miscellaneous functions

\par
Functions that don't fit into any of the other sections.


\subsubsection always

\par
A lazy identity function.
\snippet Test_TMPLDocumentation.cpp tmpl::always

\see \ref identity


\subsubsection apply

\par
Calls a \ref metalambdas "metalambda" with given arguments.
\snippet Test_TMPLDocumentation.cpp tmpl::apply


\subsubsection count

\par
Returns the number of template parameters provided as an \ref integral_constant
of `unsigned int`.
\snippet Test_TMPLDocumentation.cpp tmpl::count


\subsubsection eval_if

\par
A lazy metafunction that, if the conditional (first) argument has a true
`value`, evaluates and returns the result of the first lazy metafunction (*not*
metalambda), otherwise, evaluates and returns the result of the second lazy
metafunction.
\snippet Test_TMPLDocumentation.cpp tmpl::eval_if

\par
This is performs lazy evaluation of conditional branches outside of a
metalambda.


\subsubsection eval_if_c

\par
The same as \ref eval_if, but takes its first argument as a `bool` instead of a
type.
\snippet Test_TMPLDocumentation.cpp tmpl::eval_if_c


\subsubsection has_type

\par
A lazy metafunction that returns its second argument (defaulting to `void`),
ignoring its first argument.
\snippet Test_TMPLDocumentation.cpp tmpl::has_type

\par
This can be used to expand a parameter pack to repetitions of the same type.
\snippet Test_TMPLDocumentation.cpp tmpl::has_type:pack_expansion
\snippet Test_TMPLDocumentation.cpp tmpl::has_type:pack_expansion:asserts


\subsubsection if_

\par
A lazy metafunction that returns the second argument if the `value` static
member value of the first is true, and otherwise the third.
\snippet Test_TMPLDocumentation.cpp tmpl::if_

\warning
The second and third arguments are both evaluated, independent of which is
returned.  Use \ref defer "defer" or \ref eval_if if this is undesirable.


\subsubsection if_c

\par
The same as std::conditional.
\snippet Test_TMPLDocumentation.cpp tmpl::if_c


\subsubsection is_set

\par
Tests if all of its arguments are distinct, producing a \ref integral_constant
"bool_".
\snippet Test_TMPLDocumentation.cpp tmpl::is_set


\subsubsection repeat

\par
Calls a unary eager metafunction on a given type, then on the result of that,
then on the result of that, and so on, up to a specified (as an \ref
integral_constant) number of calls.
\snippet Test_TMPLDocumentation.cpp tmpl::repeat

\note
This function has a lazy version, but it cannot be used in a metalambda because
the template template parameter prevents manipulation of the parameter list.
\snippet Test_TMPLDocumentation.cpp tmpl::repeat:lazy


\subsubsection sizeof_

\par
A lazy metafunction that computes `sizeof` its argument as an `unsigned int`
\ref integral_constant.
\snippet Test_TMPLDocumentation.cpp tmpl::sizeof_


\subsubsection substitute

\par
Substitutes values for \ref args "args" (but *not* `_1` or `_2`) appearing in a
type.
\snippet Test_TMPLDocumentation.cpp tmpl::substitute


\subsubsection type_from

\par
Extracts the type from a \ref type_.
\snippet Test_TMPLDocumentation.cpp tmpl::type_from

\remark
This function will work on any class with a `type` type alias, but should only
be used with \ref type_ for clarity.


\subsubsection wrap

\par
Replaces the head of a sequence.
\snippet Test_TMPLDocumentation.cpp tmpl::wrap

\note
This function has a lazy version, but it cannot be used in a metalambda because
the template template parameter prevents manipulation of the parameter list.
\snippet Test_TMPLDocumentation.cpp tmpl::wrap:lazy


\subsection runtime Runtime functionality

\par
Actual C++ functions.

\par
The examples in this section use the following definition:
\snippet Test_TMPLDocumentation.cpp runtime_declarations


\subsubsection for_each_args

\par
Calls the first argument (a functor) on each of the remaining arguments, in
order.  Returns the functor.
\snippet Test_TMPLDocumentation.cpp tmpl::for_each_args:defs
\snippet Test_TMPLDocumentation.cpp tmpl::for_each_args

\note
The functor must be copyable.  This is a bug.

\par

\note
This uses a std::reference_wrapper internally, but I don't see a reason for
that.  If it were removed then this function could be constexpr (before C++20).


\subsubsection for_each

\par
Calls a functor on \ref type_ objects wrapping each type in a list, in order.
Returns the functor.
\snippet Test_TMPLDocumentation.cpp tmpl::for_each:defs
\snippet Test_TMPLDocumentation.cpp tmpl::for_each

\note
The functor must be copyable.  This is a bug.

\par

\note
An object of the list template parameter type is constructed, so the list must
be a complete type.

\see \ref type_from


\section oddities Bugs/Oddities

* \ref join has eager and lazy versions that don't agree.

* \ref push_front and \ref pop_front have lazy versions, but \ref push_back,
  and \ref pop_back do not.

* \ref reverse_range validates its arguments, but \ref range does not.
  (Probably because the former is called incorrectly more often.)

* \ref repeat and \ref wrap have unusable (in metalambdas) lazy versions
  (because they have a template template parameter).

* \ref replace has a completely broken lazy version (because it is a
  non-trivial type alias instead of a struct).

* Brigand inconsistently uses `unsigned int` and `size_t` for size-related
  things.  (Most blatantly, the result of \ref sizeof_ is represented as an
  `unsigned int`.)

* Brigand has a file containing operator overloads for \ref integral_constant
  ""s, but it is not included by the main brigand header.  They work poorly,
  mostly because it inexplicably puts them all in namespace std where the
  compiler can't find them.


\section TODO



```
    155 tmpl::conditional_t
     42 tmpl::remove_duplicates
     17 tmpl::list_difference
     14 tmpl::at
      8 tmpl::all
      8 tmpl::_state
      8 tmpl::_element
      7 tmpl::fold
      6 tmpl::make_sequence
      6 tmpl::erase
      5 tmpl::sort
      5 tmpl::found
      3 tmpl::map
      3 tmpl::has_key
      2 tmpl::get_source
      2 tmpl::get_destination
      2 tmpl::edge
      2 tmpl::any
      1 tmpl::insert
      1 tmpl::find
```

```
./adapted:
fusion.hpp        -
integral_list.hpp -
list.hpp          -
pair.hpp          -
tuple.hpp         -
variant.hpp       -

./algorithms:
all.hpp           -
any.hpp           -
count.hpp         - Done
find.hpp          -
flatten.hpp       - Done
fold.hpp          -
for_each.hpp      - Done
for_each_args.hpp - Done
index_of.hpp      - Done
is_set.hpp        - Done
merge.hpp         -
none.hpp          -
partition.hpp     -
remove.hpp        - Done
replace.hpp       - Done
reverse.hpp       - Done
select.hpp        -
sort.hpp          -
split.hpp         - Done
split_at.hpp      - Done
transform.hpp     - Done
wrap.hpp          - Done

./functions:
eval_if.hpp  - Done
if.hpp       - Done

./functions/arithmetic:
complement.hpp          - Done
divides.hpp             - Done
identity.hpp            - Done
max.hpp                 - Done
min.hpp                 - Done
minus.hpp               - Done
modulo.hpp              - Done
negate.hpp              - Done
next.hpp                - Done
plus.hpp                - Done
prev.hpp                - Done
times.hpp               - Done

./functions/bitwise:
bitand.hpp           - Done
bitor.hpp            - Done
bitxor.hpp           - Done
shift_left.hpp       - Done
shift_right.hpp      - Done

./functions/comparison:
equal_to.hpp            - Done
greater.hpp             - Done
greater_equal.hpp       - Done
less.hpp                - Done
less_equal.hpp          - Done
not_equal_to.hpp        - Done

./functions/lambda:
apply.hpp           - Done
bind.hpp            - Done
substitute.hpp      - Done

./functions/logical:
and.hpp              - Done
not.hpp              - Done
or.hpp               - Done
xor.hpp              - Done

./functions/misc:
always.hpp        - Done
repeat.hpp        - Done
sizeof.hpp        - Done

./sequences:
append.hpp             - Done
at.hpp                 -
back.hpp               - Done
clear.hpp              - Done
contains.hpp           -
erase.hpp              -
filled_list.hpp        - Done
front.hpp              - Done
has_key.hpp            -
insert.hpp             -
keys_as_sequence.hpp   -
list.hpp               - Done
make_sequence.hpp      -
map.hpp                -
pair.hpp               - Done
range.hpp              - Done
set.hpp                -
size.hpp               - Done
values_as_sequence.hpp -

./types:
args.hpp              - Done
bool.hpp              - Done
empty_base.hpp        - Used in inherit.hpp, figure out after that
has_type.hpp          - Done
inherit.hpp           -
inherit_linearly.hpp  -
integer.hpp           - Done
integral_constant.hpp - Done
no_such_type.hpp      - Done
operators.hpp         - Broken and unused
real.hpp              -
type.hpp              - Done
voidp.hpp             - Not included by main header.  Special case of has_type.
```
