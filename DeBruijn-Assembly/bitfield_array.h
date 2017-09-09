#ifndef IGEN_BASE_BITFIELD_ARRAY_H_
#define IGEN_BASE_BITFIELD_ARRAY_H_

#if defined(EXEN_DEBUG_FLAG)
  #include <cassert>
#endif

namespace IGen {
namespace base {

/// Member of a bitfield that allows indexed addressing
/// @tparam T Type of the value of the bitfield
/// @tparam BASEOFFSET Offset of the array from the starting location of the bitfield
/// @tparam BITSPERITEM Size of individual element in the array
/// @tparam NUMITEMS Number of items in the array
template <typename T, int BASEOFFSET, int BITSPERITEM, int NUMITEMS>
struct bitfield_array {
  /// The variable that stores the entire bitfield. The array corresponds to a subset of bits of this value
  T value_;

  /// Maximum value permissible in an element of the bitfield array
  static const T maximum_ = (T(1) << BITSPERITEM) - 1;

  /// Accessor function for the maximum value
  inline T maximum() const { return maximum_; }

  /// @return Number of items in the array. Returns the template parameter
  inline int num_items() const { return NUMITEMS; }

  /// Class defines each element of the array
  class element {
   private:
    /// Reference to the value of the entire bitfield
    T& value_;

    /// Offset of the element from the beginning of the bitfield
    int offset_;

   public:
    /// Construct an element by passing in the offset from the beginning of value
    /// @param value Reference to the starting of the bitfield
    /// @param offset Offset of this element from the starting of the bitfield
    inline element(T& value, int offset) : value_(value), offset_(offset) {}

    /// Mask to get only the bits corresponding to this element from the entire bitfield
    inline T mask() const { return maximum_ << offset_; }

    ///
    inline operator T() const { return (value_ >> offset_) & maximum_; }

    /// Overloaded assignment operator
    inline element& operator=(T other) {
    #if defined(EXEN_DEBUG_FLAG)
      assert(other <= maximum_);
    #endif
      value_ = (value_ & ~mask()) | (other << offset_);
      return *this;
    }

    /// Overloaded += operator
    inline element& operator+=(T other) {
    #if defined(EXEN_DEBUG_FLAG)
      assert(T(*this) + other <= maximum_);  // No overflow
    #endif
      value_ += (other << offset_);
      return *this;
    }

    /// Overloaded -= operator
    inline element& operator-=(T other) {
    #if defined(EXEN_DEBUG_FLAG)
      assert(T(*this) >= other);  // No underflow
    #endif
      value_ -= (other << offset_);
      return *this;
    }

    /// Overloaded prefix increment operator
    inline element& operator++() { return *this += 1; }

    /// Overloaded postfix increment operator
    inline element& operator++(int) { return *this += 1; }

    /// Overloaded prefix decrement operator
    inline element& operator--() { return *this -= 1; }

    /// Overloaded postfix decrement operator
    inline element& operator--(int) { return *this -= 1; }
  };

  /// Get an individual element of the array with its index
  element operator[](int i) {
#if defined(EXEN_DEBUG_FLAG)
    assert(i >= 0 && i < NUMITEMS);
#endif
    return element(value_, BASEOFFSET + BITSPERITEM*i);
  }

  /// Get an individual element of the array with its index
  const element operator[](int i) const {
#if defined(EXEN_DEBUG_FLAG)
    assert(i >= 0 && i < NUMITEMS);
#endif
    return element(value_, BASEOFFSET + BITSPERITEM*i);
  }

  // Ensure that the template arguments for this bitfield array is correct
  static_assert(BASEOFFSET + BITSPERITEM*NUMITEMS <= static_cast<int>(sizeof(T))*8,
                "Array exceeds bitfield boundaries");
  static_assert(BITSPERITEM < static_cast<int>(sizeof(T))*8, "Cannot fill entire bitfield with one element");
};

}  // namespace base
}  // namespace IGen

#endif  // IGEN_BASE_BITFIELD_ARRAY_H_
