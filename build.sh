#!/bin/bash

# Clean
make clean
# Compile
make -j8
make test -j8
# Test
bin/shocktest
