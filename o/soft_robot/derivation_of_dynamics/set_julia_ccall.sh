#!/bin/bash

echo "hoge"
mkdir -p ../derived/ikko_dake/eqs/c_so
gcc -shared -fPIC -o ../derived/ikko_dake/eqs/c_so/Phi0.so ../derived/ikko_dake/eqs/c_src/Phi0.c
gcc -shared -fPIC -o ../derived/ikko_dake/eqs/c_so/Theta0.so ../derived/ikko_dake/eqs/c_src/Theta0.c
gcc -shared -fPIC -o ../derived/ikko_dake/eqs/c_so/M.so ../derived/ikko_dake/eqs/c_src/M.c
gcc -shared -fPIC -o ../derived/ikko_dake/eqs/c_so/C.so ../derived/ikko_dake/eqs/c_src/C.c
gcc -shared -fPIC -o ../derived/ikko_dake/eqs/c_so/G.so ../derived/ikko_dake/eqs/c_src/G.c

gcc -shared -fPIC -o ../derived/ikko_dake/eqs/c_so/invM.so \
../derived/ikko_dake/eqs/c_src/invM/invM.c \
../derived/ikko_dake/eqs/c_src/invM/invM_0_0.c \
../derived/ikko_dake/eqs/c_src/invM/invM_0_1.c \
../derived/ikko_dake/eqs/c_src/invM/invM_0_2.c \
../derived/ikko_dake/eqs/c_src/invM/invM_1_0.c \
../derived/ikko_dake/eqs/c_src/invM/invM_1_1.c \
../derived/ikko_dake/eqs/c_src/invM/invM_1_2.c \
../derived/ikko_dake/eqs/c_src/invM/invM_2_0.c \
../derived/ikko_dake/eqs/c_src/invM/invM_2_1.c \
../derived/ikko_dake/eqs/c_src/invM/invM_2_2.c

gcc -shared -fPIC -o ../derived/ikko_dake/eqs/c_so/fx.so \
../derived/ikko_dake/eqs/c_src/fx/fx.c \
../derived/ikko_dake/eqs/c_src/fx/fx_0.c \
../derived/ikko_dake/eqs/c_src/fx/fx_1.c \
../derived/ikko_dake/eqs/c_src/fx/fx_2.c \
../derived/ikko_dake/eqs/c_src/fx/fx_3.c \
../derived/ikko_dake/eqs/c_src/fx/fx_4.c \
../derived/ikko_dake/eqs/c_src/fx/fx_5.c \
../derived/ikko_dake/eqs/c_src/fx/fx_6.c \
../derived/ikko_dake/eqs/c_src/fx/fx_7.c \
../derived/ikko_dake/eqs/c_src/fx/fx_8.c

gcc -shared -fPIC -o ../derived/ikko_dake/eqs/c_so/q_dot_dot.so \
../derived/ikko_dake/eqs/c_src/q_dot_dot/q_dot_dot.c \
../derived/ikko_dake/eqs/c_src/q_dot_dot/q_dot_dot_0.c \
../derived/ikko_dake/eqs/c_src/q_dot_dot/q_dot_dot_1.c \
../derived/ikko_dake/eqs/c_src/q_dot_dot/q_dot_dot_2.c \

gcc -shared -fPIC -o ../derived/ikko_dake/eqs/c_so/A.so \
../derived/ikko_dake/eqs/c_src/A/A.c \
../derived/ikko_dake/eqs/c_src/A/A_0_0.c \
../derived/ikko_dake/eqs/c_src/A/A_1_0.c \
../derived/ikko_dake/eqs/c_src/A/A_2_0.c \
../derived/ikko_dake/eqs/c_src/A/A_3_0.c \
../derived/ikko_dake/eqs/c_src/A/A_4_0.c \
../derived/ikko_dake/eqs/c_src/A/A_5_0.c \
../derived/ikko_dake/eqs/c_src/A/A_6_0.c \
../derived/ikko_dake/eqs/c_src/A/A_7_0.c \
../derived/ikko_dake/eqs/c_src/A/A_8_0.c \
../derived/ikko_dake/eqs/c_src/A/A_0_1.c \
../derived/ikko_dake/eqs/c_src/A/A_1_1.c \
../derived/ikko_dake/eqs/c_src/A/A_2_1.c \
../derived/ikko_dake/eqs/c_src/A/A_3_1.c \
../derived/ikko_dake/eqs/c_src/A/A_4_1.c \
../derived/ikko_dake/eqs/c_src/A/A_5_1.c \
../derived/ikko_dake/eqs/c_src/A/A_6_1.c \
../derived/ikko_dake/eqs/c_src/A/A_7_1.c \
../derived/ikko_dake/eqs/c_src/A/A_8_1.c \
../derived/ikko_dake/eqs/c_src/A/A_0_2.c \
../derived/ikko_dake/eqs/c_src/A/A_1_2.c \
../derived/ikko_dake/eqs/c_src/A/A_2_2.c \
../derived/ikko_dake/eqs/c_src/A/A_3_2.c \
../derived/ikko_dake/eqs/c_src/A/A_4_2.c \
../derived/ikko_dake/eqs/c_src/A/A_5_2.c \
../derived/ikko_dake/eqs/c_src/A/A_6_2.c \
../derived/ikko_dake/eqs/c_src/A/A_7_2.c \
../derived/ikko_dake/eqs/c_src/A/A_8_2.c \
../derived/ikko_dake/eqs/c_src/A/A_0_3.c \
../derived/ikko_dake/eqs/c_src/A/A_1_3.c \
../derived/ikko_dake/eqs/c_src/A/A_2_3.c \
../derived/ikko_dake/eqs/c_src/A/A_3_3.c \
../derived/ikko_dake/eqs/c_src/A/A_4_3.c \
../derived/ikko_dake/eqs/c_src/A/A_5_3.c \
../derived/ikko_dake/eqs/c_src/A/A_6_3.c \
../derived/ikko_dake/eqs/c_src/A/A_7_3.c \
../derived/ikko_dake/eqs/c_src/A/A_8_3.c \
../derived/ikko_dake/eqs/c_src/A/A_0_4.c \
../derived/ikko_dake/eqs/c_src/A/A_1_4.c \
../derived/ikko_dake/eqs/c_src/A/A_2_4.c \
../derived/ikko_dake/eqs/c_src/A/A_3_4.c \
../derived/ikko_dake/eqs/c_src/A/A_4_4.c \
../derived/ikko_dake/eqs/c_src/A/A_5_4.c \
../derived/ikko_dake/eqs/c_src/A/A_6_4.c \
../derived/ikko_dake/eqs/c_src/A/A_7_4.c \
../derived/ikko_dake/eqs/c_src/A/A_8_4.c \
../derived/ikko_dake/eqs/c_src/A/A_0_5.c \
../derived/ikko_dake/eqs/c_src/A/A_1_5.c \
../derived/ikko_dake/eqs/c_src/A/A_2_5.c \
../derived/ikko_dake/eqs/c_src/A/A_3_5.c \
../derived/ikko_dake/eqs/c_src/A/A_4_5.c \
../derived/ikko_dake/eqs/c_src/A/A_5_5.c \
../derived/ikko_dake/eqs/c_src/A/A_6_5.c \
../derived/ikko_dake/eqs/c_src/A/A_7_5.c \
../derived/ikko_dake/eqs/c_src/A/A_8_5.c \
../derived/ikko_dake/eqs/c_src/A/A_0_6.c \
../derived/ikko_dake/eqs/c_src/A/A_1_6.c \
../derived/ikko_dake/eqs/c_src/A/A_2_6.c \
../derived/ikko_dake/eqs/c_src/A/A_3_6.c \
../derived/ikko_dake/eqs/c_src/A/A_4_6.c \
../derived/ikko_dake/eqs/c_src/A/A_5_6.c \
../derived/ikko_dake/eqs/c_src/A/A_6_6.c \
../derived/ikko_dake/eqs/c_src/A/A_7_6.c \
../derived/ikko_dake/eqs/c_src/A/A_8_6.c \
../derived/ikko_dake/eqs/c_src/A/A_0_7.c \
../derived/ikko_dake/eqs/c_src/A/A_1_7.c \
../derived/ikko_dake/eqs/c_src/A/A_2_7.c \
../derived/ikko_dake/eqs/c_src/A/A_3_7.c \
../derived/ikko_dake/eqs/c_src/A/A_4_7.c \
../derived/ikko_dake/eqs/c_src/A/A_5_7.c \
../derived/ikko_dake/eqs/c_src/A/A_6_7.c \
../derived/ikko_dake/eqs/c_src/A/A_7_7.c \
../derived/ikko_dake/eqs/c_src/A/A_8_7.c \
../derived/ikko_dake/eqs/c_src/A/A_0_8.c \
../derived/ikko_dake/eqs/c_src/A/A_1_8.c \
../derived/ikko_dake/eqs/c_src/A/A_2_8.c \
../derived/ikko_dake/eqs/c_src/A/A_3_8.c \
../derived/ikko_dake/eqs/c_src/A/A_4_8.c \
../derived/ikko_dake/eqs/c_src/A/A_5_8.c \
../derived/ikko_dake/eqs/c_src/A/A_6_8.c \
../derived/ikko_dake/eqs/c_src/A/A_7_8.c \
../derived/ikko_dake/eqs/c_src/A/A_8_8.c \