#!/usr/bin/env bash

# -fPIC is position-independent code, mostly for obscure systems
# -fpermissive is worrisome

PY_VER="python2.7"
PY_LOC="/usr/include/{PY_VER}"

g++ -fPIC -fpermissive -w -c NRpyDNAcode.cpp -o NRpyDNAcode.o \
    -I"{PY_LOC}" \
    -I/usr/local/lib/"{PY_VER}"/dist-packages/numpy/core/include

g++ -shared NRpyDNAcode.o -o NRpyDNAcode.so

g++ -fPIC -fpermissive -w -c NRpyRS.cpp -o NRpyRS.o -I"{PY_LOC}" \
 -I/usr/local/lib/"{PY_VER}"/dist-packages/numpy/core/include
g++ -shared NRpyRS.o -o NRpyRS.so

# echo "done"
