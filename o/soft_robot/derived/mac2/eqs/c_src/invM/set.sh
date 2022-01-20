#!/bin/bash

gcc -shared -fPIC -o invM.so \
invM.c \
invM_0_0.c invM_0_1.c invM_0_2.c \
invM_1_0.c invM_1_1.c invM_1_2.c \
invM_2_0.c invM_2_1.c invM_2_2.c