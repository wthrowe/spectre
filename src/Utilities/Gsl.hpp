
/// \file
/// Defines functions and classes from the GSL

#pragma once

#pragma GCC system_header

// The code in this file is adapted from Microsoft's GSL that can be found at
// https://github.com/Microsoft/GSL
// The original license and copyright are:
///////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2015 Microsoft Corporation. All rights reserved.
//
// This code is licensed under the MIT License (MIT).
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//
///////////////////////////////////////////////////////////////////////////////
// The code changes are because SpECTRE is not allowed to throw under any
// circumstances that cannot be guaranteed to be caught and so all throw's
// are replaced by hard errors (ERROR).

#include "ErrorHandling/ExpectsAndEnsures.hpp"
#include "Utilities/ForceInline.hpp"
#include "Utilities/Literals.hpp"
#include "Utilities/PrettyType.hpp"

#if defined(__clang__) || defined(__GNUC__)

/*!
 * \ingroup Utilities
 * The if statement is expected to evaluate true most of the time
 */
#define LIKELY(x) __builtin_expect(!!(x), 1)

/*!
 * \ingroup Utilities
 * The if statement is expected to evaluate false most of the time
 */
#define UNLIKELY(x) __builtin_expect(!!(x), 0)

#else
/*!
 * \ingroup Utilities
 * The if statement is expected to evaluate true most of the time
 */
#define LIKELY(x) (x)

/*!
 * \ingroup Utilities
 * The if statement is expected to evaluate false most of the time
 */
#define UNLIKELY(x) (x)
#endif

/*!
 * \ingroup Utilities
 * \brief Implementations from the Guideline Support Library
 */
namespace gsl {

/*!
 * \ingroup Utilities
 * \brief Cast `u` to a type `T` where the cast may result in narrowing
 */
template <class T, class U>
SPECTRE_ALWAYS_INLINE constexpr T narrow_cast(U&& u) noexcept {
  return static_cast<T>(std::forward<U>(u));
}

namespace details {
template <class T, class U>
struct is_same_signedness
    : public std::integral_constant<
          bool, std::is_signed<T>::value == std::is_signed<U>::value> {};
}  // namespace details

/*!
 * \ingroup Utilities
 * \brief A checked version of narrow_cast() that ERRORs if the cast changed
 * the value
 */
template <class T, class U>
SPECTRE_ALWAYS_INLINE T narrow(U u) {
  T t = narrow_cast<T>(u);
  if (static_cast<U>(t) != u) {
    ERROR("Failed to cast " << u << " of type " << pretty_type::get_name<U>()
                            << " to type " << pretty_type::get_name<T>());
  }
  if (not details::is_same_signedness<T, U>::value and
      ((t < T{}) != (u < U{}))) {
    ERROR("Failed to cast " << u << " of type " << pretty_type::get_name<U>()
                            << " to type " << pretty_type::get_name<T>());
  }
  return t;
}

// @{
/*!
 * \ingroup Utilities
 * \brief Retrieve a entry from a container, with checks in Debug mode that
 * the index being retrieved is valid.
 */
template <class T, std::size_t N, typename Size>
SPECTRE_ALWAYS_INLINE constexpr T& at(std::array<T, N>& arr, Size index) {
  Expects(index >= 0 and index < narrow_cast<Size>(N));
  return arr[static_cast<std::size_t>(index)];
}

template <class Cont, typename Size>
SPECTRE_ALWAYS_INLINE constexpr const typename Cont::value_type& at(
    const Cont& cont, Size index) {
  Expects(index >= 0 and index < narrow_cast<Size>(cont.size()));
  return cont[static_cast<typename Cont::size_type>(index)];
}

template <class T, typename Size>
SPECTRE_ALWAYS_INLINE constexpr const T& at(std::initializer_list<T> cont,
                                            Size index) {
  Expects(index >= 0 and index < narrow_cast<Size>(cont.size()));
  return *(cont.begin() + index);
}
// @}
}  // namespace gsl
