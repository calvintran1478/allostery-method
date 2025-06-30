#!/bin/bash

mkdir -p ./bin

gcc -O3 -march=native -o bin/threaded_discrete_mi_max src/threaded_discrete_mi_max.c -lm -pthread
gcc -O3 -march=native -o bin/threaded_discrete_mi_avg src/threaded_discrete_mi_avg.c -lm -pthread

gcc -O3 -march=native -o bin/threaded_discrete_gcc_max src/threaded_discrete_gcc_max.c -lm -pthread
gcc -O3 -march=native -o bin/threaded_discrete_gcc_avg src/threaded_discrete_gcc_avg.c -lm -pthread

gcc -O3 -march=native -o bin/threaded_discrete_vi_min src/threaded_discrete_vi_min.c -lm -pthread
gcc -O3 -march=native -o bin/threaded_discrete_vi_avg src/threaded_discrete_vi_avg.c -lm -pthread
