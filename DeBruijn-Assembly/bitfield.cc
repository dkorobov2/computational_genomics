//#include <boost/test/unit_test.hpp>
#include <cinttypes>

#include "bitfield.h"

BEGIN_BITFIELD_TYPE(test_bitfield, uint32_t)
  ADD_BITFIELD_MEMBER(A, 0, 8)
  ADD_BITFIELD_MEMBER(B, 8, 8)
  ADD_BITFIELD_MEMBER(C, 16, 8)
  ADD_BITFIELD_MEMBER(D, 24, 8)
  ADD_BITFIELD_ARRAY(ARR, 0, 8, 4)
END_BITFIELD_TYPE()

/*BOOST_AUTO_TEST_CASE(maximum_value) {
  test_bitfield x = 0;
  BOOST_CHECK(x.A.maximum() + 1 == 0x400);
  BOOST_CHECK(x.ARR.maximum() + 1 == 0x400);
}

BOOST_AUTO_TEST_CASE(one) {
  test_bitfield x = 0;
  BOOST_CHECK(x.A.one() == 1);
}

BOOST_AUTO_TEST_CASE(mask) {
  test_bitfield x;
  BOOST_CHECK(x.ARR[0].mask() == 0x3FF);
  BOOST_CHECK(x.ARR[1].mask() == 0xFFC00);
  BOOST_CHECK(x.ARR[2].mask() == 0x3FF00000);
}

BOOST_AUTO_TEST_CASE(initialization_and_assignment) {
  test_bitfield x = 0;
  x.A = 0x2AA;
  x.B = 0x2AA;
  x.C = 0x2AA;
  BOOST_CHECK(x.A == 0x2AA);
  BOOST_CHECK(x.B == 0x2AA);
  BOOST_CHECK(x.C == 0x2AA);
  BOOST_CHECK(x.ARR[0] == 0x2AA);
  BOOST_CHECK(x.ARR[1] == 0x2AA);
  BOOST_CHECK(x.ARR[2] == 0x2AA);
  x.ARR[0] = 0x3FF;
  x.ARR[1] = 0x3FF;
  x.ARR[2] = 0x3FF;
  BOOST_CHECK(x.ARR[0] == 0x3FF);
  BOOST_CHECK(x.ARR[1] == 0x3FF);
  BOOST_CHECK(x.ARR[2] == 0x3FF);
}

BOOST_AUTO_TEST_CASE(increment) {
    test_bitfield x = 0;
    x.A = 0x1;
    x.A += 0xA;
    x.B = 0x2;
    x.B += 0xA;
    x.C = 0x3;
    x.C += 0xA;
    BOOST_CHECK(x.A == 0xB);
    BOOST_CHECK(x.B == 0xC);
    BOOST_CHECK(x.C == 0xD);
    x.ARR[0] = 0x1;
    x.ARR[0] += 0xA;
    x.ARR[1] = 0x2;
    x.ARR[1] += 0xA;
    x.ARR[2] = 0x3;
    x.ARR[2] += 0xA;
    BOOST_CHECK(x.ARR[0] == 0xB);
    BOOST_CHECK(x.ARR[1] == 0xC);
    BOOST_CHECK(x.ARR[2] == 0xD);
}

BOOST_AUTO_TEST_CASE(decrement) {
    test_bitfield x = 0;
    x.A = 0xB;
    x.A -= 0xA;
    x.B = 0xC;
    x.B -= 0xA;
    x.C = 0xD;
    x.C -= 0xA;
    BOOST_CHECK(x.A == 0x1);
    BOOST_CHECK(x.B == 0x2);
    BOOST_CHECK(x.C == 0x3);
    x.ARR[0] = 0xB;
    x.ARR[0] -= 0xA;
    x.ARR[1] = 0xC;
    x.ARR[1] -= 0xA;
    x.ARR[2] = 0xD;
    x.ARR[2] -= 0xA;
    BOOST_CHECK(x.ARR[0] == 0x1);
    BOOST_CHECK(x.ARR[1] == 0x2);
    BOOST_CHECK(x.ARR[2] == 0x3);
}*/
