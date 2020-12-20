#!/usr/bin/env bash

# -fPIC is position-independent code, mostly for obscure systems
# -fpermissive is worrisome
# TODO: create a makefile

PY_VER="python3.7"
PY_LOC="/usr/include/{PY_VER}"

warn_opts="-Wall -Werror -Wpedantic"
cpp_ver="-std=c++17"

g++ "${warn_opts}" "${cpp_ver}" -fPIC -w -c NRpyDNAcode.cpp -o NRpyDNAcode.o \
    -I"{PY_LOC}" \
    -I/usr/local/lib/"{PY_VER}"/dist-packages/numpy/core/include

g++ "${warn_opts}" "${cpp_ver}" -shared NRpyDNAcode.o -o NRpyDNAcode.so

g++ "${warn_opts}" "${cpp_ver}" -fPIC -w -c NRpyRS.cpp -o NRpyRS.o -I"{PY_LOC}" \
 -I/usr/local/lib/"{PY_VER}"/dist-packages/numpy/core/include

g++ "${warn_opts}" "${cpp_ver}" -shared NRpyRS.o -o NRpyRS.so

# echo "done"
