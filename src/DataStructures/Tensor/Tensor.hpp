// Distributed under the MIT License.
// See LICENSE.txt for details.

/// \file
/// Defines classes for Tensor

#pragma once

#include <cstddef>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Expressions/Contract.hpp"
#include "DataStructures/Tensor/Expressions/TensorExpression.hpp"
#include "DataStructures/Tensor/IndexType.hpp"
#include "DataStructures/Tensor/Structure.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "ErrorHandling/Error.hpp"
#include "Parallel/PupStlCpp11.hpp"
#include "Utilities/ForceInline.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeArray.hpp"
#include "Utilities/PrettyType.hpp"
#include "Utilities/StdHelpers.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TypeTraits.hpp"

/// \cond
template <typename X, typename Symm = Symmetry<>,
          typename IndexLs = index_list<>>
class Tensor;
/// \endcond

namespace Tensor_detail {
template <typename T, typename = void>
struct is_tensor : std::false_type {};
template <typename X, typename Symm, typename IndexLs>
struct is_tensor<Tensor<X, Symm, IndexLs>> : std::true_type {};
}  // namespace Tensor_detail

/*!
 * \ingroup Tensor
 * \brief Represents an object with multiple components
 *
 * \details
 * Tensor is a container that represents indexable geometric objects. Each index
 * has its own dimension, valence, and frame and must be either spatial or
 * spacetime. Note that the dimension passed to `SpatialIndex` and
 * `SpacetimeIndex` is always the spatial dimension of the index. Tensors with
 * symmetric indices are stored only once and must be of the same
 * type. A list of available type aliases can be found in the ::tnsr namespace
 * where the adopted conventions are:
 *
 * 1. Upper case for contravariant or upper indices, lower case for covariant or
 * lower indices.
 *
 * 2. `a, b, c, d` are used for spacetime indices while `i, j, k, l` are used
 * for spatial indices.
 *
 * 3. `::Scalar` is not inside the `::tnsr` namespace but is used to represent
 * a scalar with no indices.
 *
 * \example
 * \snippet Test_Tensor.cpp scalar
 * \snippet Test_Tensor.cpp spatial_vector
 * \snippet Test_Tensor.cpp spacetime_vector
 * \snippet Test_Tensor.cpp rank_3_122
 *
 * \tparam X the type held
 * \tparam Symm the ::Symmetry of the indices
 * \tparam IndexLs a typelist of \ref SpacetimeIndex "TensorIndexType"'s
 */
template <typename X, typename Symm, template <typename...> class IndexLs,
          typename... Indices>
class Tensor<X, Symm, IndexLs<Indices...>> {
  /// The number of \ref SpacetimeIndex "TensorIndexType"'s
  ///
  /// Note: Scalars need to have 1 so we can still store their data.
  static constexpr auto num_tensor_indices =
      sizeof...(Indices) == 0 ? 1 : sizeof...(Indices);

  /// The Tensor_detail::Structure for the particular tensor index structure
  ///
  /// Each tensor index structure, e.g. \f$T_{ab}\f$, \f$T_a{}^b\f$ or
  /// \f$T^{ab}\f$ has its own Tensor_detail::TensorStructure that holds
  /// information about how the data is stored, what the multiplicity of the
  /// stored indices are, the number of (independent) components, etc.
  using structure = Tensor_detail::Structure<Symm, Indices...>;

 public:
  /// The type of the sequence that holds the data
  using storage_type =
      std::array<X, Tensor_detail::Structure<Symm, Indices...>::size()>;
  /// The type that is stored by the Tensor
  using type = X;
  /// Typelist of the symmetry of the Tensor
  ///
  /// \details
  /// For a rank-3 tensor symmetric in the last two indices,
  /// \f$T_{a(bc)}\f$, the ::Symmetry is `<2, 1, 1>`. For a non-symmetric rank-2
  /// tensor the ::Symmetry is `<2, 1>`.
  using symmetry = Symm;
  /// Typelist of the \ref SpacetimeIndex "TensorIndexType"'s that the
  /// Tensor has
  using index_list = typelist<Indices...>;
  /// The type of the TensorExpression that would represent this Tensor in a
  /// tensor expression.
  template <typename ArgsLs>
  using TE = TensorExpression<Tensor<X, Symm, tmpl::list<Indices...>>, X, Symm,
                              tmpl::list<Indices...>, ArgsLs>;

  Tensor() = default;
  ~Tensor() = default;
  Tensor(const Tensor&) = default;
  Tensor(Tensor&&) noexcept = default;  // NOLINT
  Tensor& operator=(const Tensor&) = default;
  Tensor& operator=(Tensor&&) noexcept = default;  // NOLINT

  /// \cond HIDDEN_SYMBOLS
  /// Constructor from a TensorExpression.
  ///
  /// \tparam LhsIndices the indices on the LHS of the tensor expression
  /// \tparam T the type of the TensorExpression
  /// \param tensor_expression the tensor expression being evaluated
  template <typename... LhsIndices, typename T,
            std::enable_if_t<std::is_base_of<Expression, T>::value>* = nullptr>
  Tensor(const T& tensor_expression, tmpl::list<LhsIndices...> /*meta*/) {
    static_assert(
        sizeof...(LhsIndices) == sizeof...(Indices),
        "When calling evaluate<...>(...) you must pass the same "
        "number of indices as template parameters as there are free "
        "indices on the resulting tensor. For example, auto F = "
        "evaluate<_a_t, _b_t>(G); if G has 2 free indices and you want "
        "the LHS of the equation to be F_{ab} rather than F_{ba}.");
    for (size_t i = 0; i < size(); ++i) {
      data_[i] = tensor_expression.template get<LhsIndices...>(  // NOLINT
          structure::get_canonical_tensor_index(i));
    }
  }
  /// \endcond

  /// Initialize a Vector with the value in the std::initializer_list
  ///
  /// \example
  /// \snippet Test_Tensor.cpp init_vector
  /// \param data the values of the individual components of the Vector
  template <size_t NumberOfIndices = sizeof...(Indices),
            std::enable_if_t<(NumberOfIndices == 1)>* = nullptr>
  explicit Tensor(storage_type data);

  /// Constructor that passes "args" to constructor of X and initializes each
  /// component to be the same
  template <
      typename... Args,
      std::enable_if_t<not(
          cpp17::disjunction_v<std::is_same<
              Tensor<X, Symm, IndexLs<Indices...>>, std::decay_t<Args>>...> and
          sizeof...(Args) == 1)>* = nullptr>
  explicit Tensor(Args&&... args);

  using value_type = typename storage_type::value_type;
  using reference = typename storage_type::reference;
  using const_reference = typename storage_type::const_reference;
  using iterator = typename storage_type::iterator;
  using const_iterator = typename storage_type::const_iterator;
  using pointer = typename storage_type::pointer;
  using const_pointer = typename storage_type::const_pointer;
  using reverse_iterator = typename storage_type::reverse_iterator;
  using const_reverse_iterator = typename storage_type::const_reverse_iterator;

  iterator begin() { return data_.begin(); }
  const_iterator begin() const { return data_.begin(); }
  const_iterator cbegin() const { return data_.begin(); }

  iterator end() { return data_.end(); }
  const_iterator end() const { return data_.end(); }
  const_iterator cend() const { return data_.end(); }

  reverse_iterator rbegin() { return data_.rbegin(); }
  const_reverse_iterator rbegin() const { return data_.rbegin(); }
  const_reverse_iterator crbegin() const { return data_.rbegin(); }

  reverse_iterator rend() { return data_.rend(); }
  const_reverse_iterator rend() const { return data_.rend(); }
  const_reverse_iterator crend() const { return data_.rend(); }

  // @{
  /// Get data entry using an array representing a tensor index
  ///
  /// \details
  /// Let \f$T_{abc}\f$ be a Tensor.
  /// Then `get({{0, 2, 1}})` returns the \f$T_{0 2 1}\f$ component.
  /// \param tensor_index the index at which to get the data
  template <typename T>
  SPECTRE_ALWAYS_INLINE constexpr reference get(
      const std::array<T, num_tensor_indices>& tensor_index) {
    return gsl::at(data_, structure::get_storage_index(tensor_index));
  }
  template <typename T>
  SPECTRE_ALWAYS_INLINE constexpr const_reference get(
      const std::array<T, num_tensor_indices>& tensor_index) const {
    return gsl::at(data_, structure::get_storage_index(tensor_index));
  }
  // @}
  // @{
  /// Get data entry using a list of integers representing a tensor index
  ///
  /// \details
  /// Let \f$T_{abc}\f$ be a Tensor.
  /// Then `get(0, 2, 1)` returns the \f$T_{0 2 1}\f$ component.
  /// \param n the index at which to get the data
  template <typename... N>
  constexpr reference get(N... n) {
    static_assert(
        sizeof...(Indices) == sizeof...(N),
        "the number of tensor indices specified must match the rank of "
        "the tensor");
    return gsl::at(data_, structure::get_storage_index(n...));
  }
  template <typename... N>
  constexpr const_reference get(N... n) const {
    static_assert(
        sizeof...(Indices) == sizeof...(N),
        "the number of tensor indices specified must match the rank of "
        "the tensor");
    return gsl::at(data_, structure::get_storage_index(n...));
  }
  // @}

  // @{
  /// Retrieve the index `N...` by computing the storage index at compile time
  template <int... N>
  SPECTRE_ALWAYS_INLINE constexpr reference get() {
    static_assert(sizeof...(Indices) == sizeof...(N),
                  "the number of tensor indices specified must match the rank "
                  "of the tensor");
    return gsl::at(data_, structure::template get_storage_index<N...>());
  }
  template <int... N>
  SPECTRE_ALWAYS_INLINE constexpr const_reference get() const {
    static_assert(sizeof...(Indices) == sizeof...(N),
                  "the number of tensor indices specified must match the rank "
                  "of the tensor");
    return gsl::at(data_, structure::template get_storage_index<N...>());
  }
  // @}

  // @{
  /// Retrieve a TensorExpression object with the index structure passed in
  template <typename... N,
            std::enable_if_t<cpp17::conjunction_v<tt::is_tensor_index<N>...> and
                             tmpl::is_set<N...>::value>* = nullptr>
  SPECTRE_ALWAYS_INLINE constexpr TE<tmpl::list<N...>> operator()(
      N... /*meta*/) const {
    return TE<tmpl::list<N...>>(*this);
  }

  template <typename... N,
            std::enable_if_t<cpp17::conjunction_v<tt::is_tensor_index<N>...> and
                             not tmpl::is_set<N...>::value>* = nullptr>
  SPECTRE_ALWAYS_INLINE constexpr auto operator()(N... /*meta*/) const
      -> decltype(
          TensorExpressions::fully_contracted<
              TE, replace_indices<tmpl::list<N...>, repeated<tmpl::list<N...>>>,
              tmpl::int32_t<0>,
              tmpl::size<repeated<tmpl::list<N...>>>>::apply(*this)) {
    using args_list = tmpl::list<N...>;
    using repeated_args_list = repeated<args_list>;
    return TensorExpressions::fully_contracted<
        TE, replace_indices<args_list, repeated_args_list>, tmpl::int32_t<0>,
        tmpl::size<repeated_args_list>>::apply(*this);
  }
  // @}

  // @{
  /// Return i'th component of storage vector
  constexpr reference operator[](const size_t storage_index) {
    return gsl::at(data_, storage_index);
  }
  constexpr const_reference operator[](const size_t storage_index) const {
    return gsl::at(data_, storage_index);
  }
  // @}

  /// Return the number of independent components of the Tensor
  ///
  /// \details
  /// Returns the number of independent components of the Tensor taking into
  /// account symmetries. For example, let \f$T_{ab}\f$ be a n-dimensional
  /// rank-2 symmetric tensor, then the number of independent components is
  /// \f$n(n+1)/2\f$.
  SPECTRE_ALWAYS_INLINE static constexpr size_t size() noexcept {
    return structure::size();
  }

  /// Returns the rank of the Tensor
  ///
  /// \details
  /// The rank of a tensor is the number of indices it has. For example, the
  /// tensor \f$v^a\f$ is rank-1, the tensor \f$\phi\f$ is rank-0, and the
  /// tensor \f$T_{abc}\f$ is rank-3.
  SPECTRE_ALWAYS_INLINE static constexpr size_t rank() noexcept {
    return sizeof...(Indices);
  }

  // @{
  /// Given an iterator or storage index, get the canonical tensor index.
  /// For scalars this is defined to be std::array<int, 1>{{0}}
  SPECTRE_ALWAYS_INLINE constexpr std::array<
      size_t, sizeof...(Indices) == 0 ? 1 : sizeof...(Indices)>
  get_tensor_index(const const_iterator& iter) noexcept {
    return structure::get_canonical_tensor_index(
        static_cast<size_t>(iter - begin()));
  }
  SPECTRE_ALWAYS_INLINE static constexpr std::array<
      size_t, sizeof...(Indices) == 0 ? 1 : sizeof...(Indices)>
  get_tensor_index(const size_t storage_index) noexcept {
    return structure::get_canonical_tensor_index(storage_index);
  }
  // @}

  // @{
  /// Given an iterator or storage index, get the multiplicity of an index
  ///
  /// \see TensorMetafunctions::compute_multiplicity
  SPECTRE_ALWAYS_INLINE constexpr size_t multiplicity(
      const iterator& iter) const noexcept {
    return structure::multiplicity(static_cast<size_t>(iter - begin()));
  }
  SPECTRE_ALWAYS_INLINE static constexpr size_t multiplicity(
      const size_t storage_index) noexcept {
    return structure::multiplicity(storage_index);
  }
  // @}

  // @{
  /// Get dimensionality of i'th tensor index
  ///
  /// \details
  /// Consider the tensor \f$T_{abc}\f$ where the dimensionality of each index
  /// is \f$(3, 2, 4)\f$, respectively. Then
  /// \code
  /// dimensionality(0) = 3;
  /// dimensionality(1) = 2;
  /// dimensionality(2) = 4;
  /// \endcode
  SPECTRE_ALWAYS_INLINE static constexpr size_t index_dim(
      const size_t i) noexcept {
    static_assert(sizeof...(Indices),
                  "A scalar does not have any indices from which you can "
                  "retrieve the dimensionality.");
    return gsl::at(structure::dims(), i);
  }
  template <int I>
  SPECTRE_ALWAYS_INLINE static constexpr size_t index_dim() noexcept {
    return structure::template dim<I>();
  }
  // @}

  //@{
  /// Return an array corresponding to the ::Symmetry of the Tensor
  SPECTRE_ALWAYS_INLINE static constexpr std::array<int, sizeof...(Indices)>
  symmetries() noexcept {
    return structure::symmetries();
  }
  //@}

  //@{
  /// Return array of the ::IndexType's (spatial or spacetime)
  SPECTRE_ALWAYS_INLINE static constexpr std::array<IndexType,
                                                    sizeof...(Indices)>
  index_types() noexcept {
    return structure::index_types();
  }
  //@}

  //@{
  /// Return array of dimensionality of each index
  SPECTRE_ALWAYS_INLINE static constexpr std::array<size_t, sizeof...(Indices)>
  index_dims() noexcept {
    return structure::dims();
  }
  //@}

  //@{
  /// Return array of the valence of each index (::UpLo)
  SPECTRE_ALWAYS_INLINE static constexpr std::array<UpLo, sizeof...(Indices)>
  index_valences() noexcept {
    return structure::index_valences();
  }
  //@}

  //@{
  /// Returns std::tuple of the ::Frame of each index
  SPECTRE_ALWAYS_INLINE static constexpr auto index_frames() noexcept {
    return Tensor_detail::Structure<Symm, Indices...>::index_frames();
  }
  //@}

  //@{
  /// Given a tensor index, get the canonical label associated with the
  /// canonical \ref SpacetimeIndex "TensorIndexType"
  template <typename T = int>
  static std::string component_name(
      const std::array<T, rank()>& tensor_index = std::array<T, rank()>{},
      const std::array<std::string, rank()>& axis_labels =
          make_array<rank()>(std::string(""))) {
    return structure::component_name(tensor_index, axis_labels);
  }
  //@}

  /// Copy tensor data into an `std::vector<X>` along with the
  /// component names into a `std::vector<std::string>`
  /// \requires `std::is_same<X, DataVector>::%value` is true
  std::pair<std::vector<std::string>, std::vector<X>> get_vector_of_data()
      const;

  /// \cond HIDDEN_SYMBOLS
  /// Serialization function used by Charm++
  void pup(PUP::er& p) {  // NOLINT
    p | data_;
  }
  /// \endcond

 private:
  storage_type data_;
};

// ================================================================
// Template Definitions - Variadic templates must be in header
// ================================================================

template <typename X, typename Symm, template <typename...> class IndexLs,
          typename... Indices>
template <size_t NumberOfIndices, std::enable_if_t<(NumberOfIndices == 1)>*>
Tensor<X, Symm, IndexLs<Indices...>>::Tensor(storage_type data)
    : data_(std::move(data)) {}

// The cpp17::disjunction is used to prevent the compiler from matching this
// function when it should select the move constructor.
template <typename X, typename Symm, template <typename...> class IndexLs,
          typename... Indices>
template <
    typename... Args,
    std::enable_if_t<not(
        cpp17::disjunction_v<std::is_same<Tensor<X, Symm, IndexLs<Indices...>>,
                                          std::decay_t<Args>>...> and
        sizeof...(Args) == 1)>*>
Tensor<X, Symm, IndexLs<Indices...>>::Tensor(Args&&... args)
    : data_(make_array<size()>(X(std::forward<Args>(args)...))) {}

template <typename X, typename Symm, template <typename...> class IndexLs,
          typename... Indices>
std::pair<std::vector<std::string>, std::vector<X>>
Tensor<X, Symm, IndexLs<Indices...>>::get_vector_of_data() const {
  std::vector<value_type> serialized_tensor(size());
  std::vector<std::string> component_names(size());
  for (size_t i = 0; i < data_.size(); ++i) {
    component_names[i] = component_name(get_tensor_index(i));
    serialized_tensor[i] = data_[i];
  }
  return std::make_pair(component_names, serialized_tensor);
}

template <typename X, typename Symm, template <typename...> class IndexLs,
          typename... Indices>
bool operator==(const Tensor<X, Symm, IndexLs<Indices...>>& lhs,
                const Tensor<X, Symm, IndexLs<Indices...>>& rhs) {
  return std::equal(lhs.begin(), lhs.end(), rhs.begin());
}
template <typename X, typename Symm, template <typename...> class IndexLs,
          typename... Indices>
bool operator!=(const Tensor<X, Symm, IndexLs<Indices...>>& lhs,
                const Tensor<X, Symm, IndexLs<Indices...>>& rhs) {
  return not(lhs == rhs);
}

// We place the stream operators in the header file so they do not need to be
// explicitly instantiated.
template <typename X, typename Symm, template <typename...> class IndexLs,
          typename... Indices>
std::ostream& operator<<(std::ostream& os,
                         const Tensor<X, Symm, IndexLs<Indices...>>& x) {
  static_assert(tt::is_streamable_v<decltype(os), X>,
                "operator<< is not defined for the type you are trying to "
                "stream in Tensor");
  os << "--Symmetry:  " << x.symmetries() << "\n";
  os << "--Types:     " << x.index_types() << "\n";
  os << "--Dims:      " << x.index_dims() << "\n";
  os << "--Locations: " << x.index_valences() << "\n";
  os << "--Frames:    " << x.index_frames() << "\n";
  for (size_t i = 0; i < x.size() - 1; ++i) {
    os << " T" << x.get_tensor_index(i) << "=" << x[i]
       << "\n     Multiplicity: " << x.multiplicity(i) << " Index: " << i
       << "\n";
  }
  size_t i = x.size() - 1;
  os << " T" << x.get_tensor_index(i) << "=" << x[i]
     << "\n     Multiplicity: " << x.multiplicity(i) << " Index: " << i;
  return os;
}
