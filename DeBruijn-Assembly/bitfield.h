#ifndef IGEN_BASE_BITFIELD_H_
#define IGEN_BASE_BITFIELD_H_

#include "bitfield_member.h"
#include "bitfield_array.h"

#define BEGIN_BITFIELD_TYPE(NAME, T) \
    union NAME \
    { \
        struct wrapper { T value_; } wrapper_; \
        NAME(T val = 0) { wrapper_.value_ = val; } \
        NAME& operator=(T val) { wrapper_.value_ = val; return *this; } \
        operator T&() { return wrapper_.value_; } \
        operator T() const { return wrapper_.value_; } \
        using storage_type = T;

#define ADD_BITFIELD_MEMBER(memberName, offset, bits) \
        IGen::base::bitfield_member<storage_type, offset, bits> memberName;

#define ADD_BITFIELD_ARRAY(memberName, offset, bits, numItems) \
        IGen::base::bitfield_array<storage_type, offset, bits, numItems> memberName;

#define END_BITFIELD_TYPE() \
    };

#endif  // IGEN_BASE_BITFIELD_H_
