#!/bin/bash

set -ex

export GK_SYSTEM=gnu_ubuntu
make -I Makefiles -j2

make -I Makefiles build-pfunit-library
make -I Makefiles check
