#!/bin/bash

mkdir -p ./derived/so

gcc -shared -fPIC -o ./derived/so/Phi_0.so ./derived/c_src/Phi_s/Phi_0.c
gcc -shared -fPIC -o ./derived/so/Phi_1.so ./derived/c_src/Phi_s/Phi_1.c
gcc -shared -fPIC -o ./derived/so/Phi_2.so ./derived/c_src/Phi_s/Phi_2.c
gcc -shared -fPIC -o ./derived/so/Phi_3.so ./derived/c_src/Phi_s/Phi_3.c
gcc -shared -fPIC -o ./derived/so/Phi_4.so ./derived/c_src/Phi_s/Phi_4.c