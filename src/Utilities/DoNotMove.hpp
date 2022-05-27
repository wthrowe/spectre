// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

/// \ingroup UtilitiesGroup
/// The "opposite" of std::move.  Casts an rvalue (or anything else)
/// to a const lvalue reference, which will prevent its use as the
/// argument to a move constructor or assignment operator.
template <typename T>
const T& do_not_move(T&& t) { return t; }
