#!/bin/bash

set -ex

cmake . -B build -DSTELLA_ENABLE_TESTS=on
cmake --build build -j2 --target check
