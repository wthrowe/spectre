// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "DataStructures/DataVector.hpp"

#include <pup_stl.h>

#include "Utilities/StdHelpers.hpp"

DataVector::DataVector(const size_t size, const double value)
    : size_(size),
      owned_data_(size_, value),
      data_(owned_data_.data(), size_) {}

DataVector::DataVector(std::initializer_list<double> list)
    : size_(list.size()), owned_data_(list), data_(owned_data_.data(), size_) {}

DataVector::DataVector(double* start, size_t size)
    : size_(size), owned_data_(0), data_(start, size_), owning_(false) {}

/// \cond HIDDEN_SYMBOLS
DataVector::DataVector(const DataVector& rhs) {
  size_ = rhs.size();
  if (rhs.is_owning()) {
    owned_data_ = rhs.owned_data_;
  } else {
    owned_data_ = InternalStorage_t(rhs.begin(), rhs.end());
  }
  data_ = decltype(data_){owned_data_.data(), size_};
}

DataVector& DataVector::operator=(const DataVector& rhs) {
  if (this == &rhs) {
    return *this;
  }
  if (owning_) {
    size_ = rhs.size();
    if (rhs.is_owning()) {
      owned_data_ = rhs.owned_data_;
    } else {
      owned_data_ = InternalStorage_t(rhs.begin(), rhs.end());
    }
    data_ = decltype(data_){owned_data_.data(), size_};
  } else {
    ASSERT(rhs.size() == size(), "Must copy into same size");
    std::copy(rhs.begin(), rhs.end(), begin());
  }
  return *this;
}

DataVector& DataVector::operator=(DataVector&& rhs) noexcept {
  if (this == &rhs) {
    return *this;
  }
  if (owning_ or size_ == 0) {
    size_ = rhs.size_;
    owned_data_ = std::move(rhs.owned_data_);
    data_ = std::move(rhs.data_);
    owning_ = rhs.owning_;
  } else {
    ASSERT(rhs.size() == size(), "Must copy into same size");
    std::copy(rhs.begin(), rhs.end(), begin());
  }
  return *this;
}
/// \endcond

void DataVector::set_data_ref(double* start, size_t size) {
  size_ = size;
  owned_data_ = decltype(owned_data_){};
  data_ = decltype(data_){start, size_};
  owning_ = false;
}

void DataVector::pup(PUP::er& p) {  // NOLINT
  p | size_;
  if (p.isUnpacking()) {
    owning_ = true;
    p | owned_data_;
    data_ = decltype(data_){owned_data_.data(), size_};
  } else {
    if (not owning_) {
      owned_data_ =
          InternalStorage_t(data_.data(), data_.data() + size_);  // NOLINT
      p | owned_data_;
      owned_data_.clear();
    } else {
      p | owned_data_;
    }
  }
}

std::ostream& operator<<(std::ostream& os, const DataVector& d) {
  // This function is inside the detail namespace StdHelpers.hpp
  StdHelpers_detail::print_helper(os, d.begin(), d.end());
  return os;
}

/// Equivalence operator for DataVector
bool operator==(const DataVector& lhs, const DataVector& rhs) {
  return lhs.size() == rhs.size() and
         std::equal(lhs.begin(), lhs.end(), rhs.begin());
}

/// Inequivalence operator for DataVector
bool operator!=(const DataVector& lhs, const DataVector& rhs) {
  return not(lhs == rhs);
}
