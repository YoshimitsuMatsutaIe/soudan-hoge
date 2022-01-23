#!/bin/bash

gcc -shared -fPIC -o q_dot_dot.so \
q_dot_dot.c \
q_dot_dot_0.c q_dot_dot_1.c q_dot_dot_2.c 