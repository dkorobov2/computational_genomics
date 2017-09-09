#ifndef IGEN_BASE_BITFIELD_MEMBER_H_
#define IGEN_BASE_BITFIELD_MEMBER_H_

#if defined(EXEN_DEBUG_FLAG)
  #include <cassert>
#endif

namespace IGen {
namespace base {

/// An individual member of a bitfield
/// @tparam T Type for the value of bitfield
/// @tparam OFFSET Base Offset inside the bitfield
/// @tparam BITS
template <typename T, int OFFSET, int BITS>
struct bitfield_member {
  /// Variable holding the bitfield. The template arguments are used to compute on a subset of bits
  T value_;

  /// Maximum value of the bitfield member
  static const T maximum_ = (T(1) << BITS) - 1;

  /// Mask to get the only the value of the bitfield member from the bitfield
  static const T mask_ = maximum_ << OFFSET;

  /// Accessor function for the maximum value of the bitfield_member
  inline T maximum() const { return maximum_; }

  /// Accessor function for the unit value for the bitfield member
  inline T one() const { return T(1) << OFFSET; }

  /// Overloaded = operator
  inline bitfield_member& operator=(T other) {
#if defined(EXEN_DEBUG_FLAG)
    assert(other <= maximum_);  // Make sure that other fits inside the bitfield
#endif
    value_ = (value_ & ~mask_) | (other << OFFSET);
    return *this;
  }

  /// Overloaded == operator
  inline bool operator==(T other) {
#if defined(EXEN_DEBUG_FLAG)
    assert(other <= maximum_);  // Make sure that other fits inside the bitfield
#endif
    return (value_ & mask_) == (other << OFFSET);
  }

  /// Overloaded += operator
  inline bitfield_member& operator+=(T other) {
#if defined(EXEN_DEBUG_FLAG)
    assert(((value_ & mask_) >> OFFSET) + other <= maximum_);  // No overflow
#endif
    value_ += (other << OFFSET);
    return *this;
  }

  /// Overloaded -= operator
  inline bitfield_member& operator-=(T other) {
#if defined(EXEN_DEBUG_FLAG)
    assert(((value_ & mask_) >> OFFSET) >= other);  // No underflow
#endif
    value_ -= (other << OFFSET);
    return *this;
  }

  /// Overloaded ++ prefix operator for the bitfield using the += operator
  inline bitfield_member& operator++() {
    return *this += 1;
  }

  /// Overloaded ++ postfix operator for the bitfield using the += operator
  inline bitfield_member& operator++(int) {
    return *this += 1;
  }

  /// Overloaded -- prefix operator for the bitfield using the -= operator
  inline bitfield_member& operator--() {
    return *this -= 1;
  }

  /// Overloaded -- postfix operator for the bitfield using the -= operator
  inline bitfield_member& operator--(int) {
    return *this -= 1;
  }

  // Ensure that the template arguments for this bitfield member is correct
  static_assert(OFFSET + BITS <= static_cast<int>(sizeof(T))*8, "Member exceeds bitfield boundaries");
  static_assert(BITS < static_cast<int>(sizeof(T))*8, "Cannot fill entire bitfield with one member");
};

}  // namespace base
}  // namespace IGen

#endif  // IGEN_BASE_BITFIELD_MEMBER_H_
