#!/bin/bash

mkdir -p ./derived/so

gcc -shared -fPIC -o ./derived/so/J_0.so ./derived/c_src/J_s/J_0.c
gcc -shared -fPIC -o ./derived/so/J_1.so ./derived/c_src/J_s/J_1.c
gcc -shared -fPIC -o ./derived/so/J_2.so ./derived/c_src/J_s/J_2.c
gcc -shared -fPIC -o ./derived/so/J_3.so ./derived/c_src/J_s/J_3.c
gcc -shared -fPIC -o ./derived/so/J_4.so ./derived/c_src/J_s/J_4.c