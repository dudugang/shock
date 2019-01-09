#!/bin/bash

g++ -std=c++11 -I lib/eigen Algebra.cpp Algebra.h Flux.cpp Flux.h Flow.cpp Flow.h main.cpp -o flow
