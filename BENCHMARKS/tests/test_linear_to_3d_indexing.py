#!/usr/bin/env python3
"""Regression tests for converting a linear index to three dimensions."""

import math
import struct
import unittest


FLOAT32_LIMIT = 2**24


def decode_integer(index, nr, nphi):
    """Decode a one-based linear index using integer arithmetic."""
    i = (index - 1) % nr + 1
    j = (index - 1) % (nr * nphi) // nr + 1
    k = (index - 1) // (nr * nphi) + 1
    return i, j, k


def encode(i, j, k, nr, nphi):
    """Encode one-based three-dimensional indices as a linear index."""
    return i + (j - 1) * nr + (k - 1) * nr * nphi


def float32(value):
    """Round a Python number to the default Fortran REAL representation."""
    return struct.unpack("f", struct.pack("f", value))[0]


def decode_k_with_real(index, nr, nphi):
    """Reproduce the former CEILING(REAL(...)/REAL(...)) expression."""
    quotient = float32(float32(index) / float32(nr * nphi))
    return math.ceil(quotient)


class TestLinearTo3DIndexing(unittest.TestCase):
    def test_integer_decode_round_trips_across_float32_limit(self):
        nr = 384
        nphi = 192
        first = FLOAT32_LIMIT - 4096
        last = FLOAT32_LIMIT + 65536

        for index in range(first, last + 1):
            decoded = decode_integer(index, nr, nphi)
            self.assertEqual(encode(*decoded, nr, nphi), index)

    def test_previous_real_decode_loses_an_index(self):
        nr = 384
        nphi = 192
        index = 16809985

        integer_k = decode_integer(index, nr, nphi)[2]
        self.assertEqual(integer_k, 229)
        self.assertEqual(decode_k_with_real(index, nr, nphi), 228)


if __name__ == "__main__":
    unittest.main()
