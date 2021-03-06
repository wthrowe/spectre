// Distributed under the MIT License.
// See LICENSE.txt for details.

/// \file
/// Defines Parallel::printf for writing to stdout

#pragma once

#include <charm++.h>
#include <sstream>
#include <string>
#include <type_traits>

#include "Utilities/TypeTraits.hpp"

namespace Parallel {
namespace detail {
/*!
 * Fundamentals and pointers are already printable so there is no conversion
 * to a std::string necessary.
 */
template <typename T,
          typename std::enable_if_t<
              std::is_fundamental<std::decay_t<
                  std::remove_pointer_t<std::decay_t<T>>>>::value or
              std::is_pointer<T>::value or
              std::is_pointer<std::decay_t<T>>::value>* = nullptr>
inline constexpr T stream_object_to_string(T&& t) {
  return t;
}

/*!
 * Stream an object of type `T` into a std::stringstream and return it as a
 * std::string so that it is printable by calling `.c_str()` on it.
 * We need a 2-phase call so that the std::string doesn't go out of scope before
 * the C-style string is passed to printf.
 */
template <
    typename T,
    typename std::enable_if_t<std::is_class<std::decay_t<T>>::value>* = nullptr>
inline std::string stream_object_to_string(T&& t) {
  static_assert(tt::is_streamable<std::stringstream, T>::value,
                "Cannot stream type and therefore it cannot be printed. Please "
                "define a stream operator for the type.");
  std::stringstream ss;
  ss << t;
  return ss.str();
}

/*!
 * Fundamentals are already printable, so nothing to do.
 */
template <typename T,
          typename std::enable_if_t<std::is_fundamental<std::decay_t<
              std::remove_pointer_t<std::decay_t<T>>>>::value>* = nullptr>
inline constexpr T get_printable_type(T&& t) {
  return t;
}

/*!
 * Get the pointer of the std::string so it can be passed to CkPrintf which
 * only works on fundamentals
 */
inline const typename std::string::value_type* get_printable_type(
    const std::string& t) {
  return t.c_str();
}

template <typename... Ts>
inline void print_helper(const std::string& format, Ts&&... t) {
  CkPrintf(format.c_str(),                               // NOLINT
           get_printable_type(std::forward<Ts>(t))...);  // NOLINT
}
}  // namespace detail

/*!
   * \ingroup Parallel
   * \brief Print an atomic message to stdout with C printf usage.
   *
   * Similar to Python, you can print any object that's streamable by passing it
   * in as an argument and using the formatter "%s". For example,
   * \code
   * std::vector<double> a{0.8, 73, 9.8};
   * Parallel::printf("%s\n", a);
   * \endcode
   */
template <typename... Args>
inline void printf(const std::string& format, Args&&... args) {
  detail::print_helper(
      format, detail::stream_object_to_string(std::forward<Args>(args))...);
}
}  // namespace Parallel
