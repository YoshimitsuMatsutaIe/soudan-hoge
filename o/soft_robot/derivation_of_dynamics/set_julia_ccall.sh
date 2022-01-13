#!/bin/bash

echo "hoge"
mkdir -p ../derived/ikko_dake/eqs/c_so
gcc -shared -fPIC -o ../derived/ikko_dake/eqs/c_so/Phi0.so ../derived/ikko_dake/eqs/c_src/Phi0.c
gcc -shared -fPIC -o ../derived/ikko_dake/eqs/c_so/Theta0.so ../derived/ikko_dake/eqs/c_src/Theta0.c
gcc -shared -fPIC -o ../derived/ikko_dake/eqs/c_so/M.so ../derived/ikko_dake/eqs/c_src/M.c
gcc -shared -fPIC -o ../derived/ikko_dake/eqs/c_so/C.so ../derived/ikko_dake/eqs/c_src/C.c
gcc -shared -fPIC -o ../derived/ikko_dake/eqs/c_so/G.so ../derived/ikko_dake/eqs/c_src/G.c
