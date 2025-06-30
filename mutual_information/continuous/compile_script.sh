#!/bin/bash

mkdir -p ./bin

gcc -O3 -march=native -o bin/threaded_continuous_gcc_avg src/threaded_continuous_gcc_avg.c -lm -pthread
