#include "fx.h"


#include "fx_0.h"
#include "fx_1.h"
#include "fx_2.h"
#include "fx_3.h"
#include "fx_4.h"
#include "fx_5.h"
#include "fx_6.h"
#include "fx_7.h"
#include "fx_8.h"

void fx(double h1, double h2, double h3, double l1, double l1_dot, double l2, double l2_dot, double l3, double l3_dot, double *out) {
    out[0] = fx_0(l1_dot);
    out[1] = fx_1(l2_dot);
    out[2] = fx_2(l3_dot);
    out[3] = fx_3(h1, h2, h3, l1, l1_dot, l2, l2_dot, l3, l3_dot);
    out[4] = fx_4(h1, h2, h3, l1, l1_dot, l2, l2_dot, l3, l3_dot);
    out[5] = fx_5(h1, h2, h3, l1, l1_dot, l2, l2_dot, l3, l3_dot);
    out[6] = fx_6(l1);
    out[7] = fx_7(l2);
    out[8] = fx_8(l3);
}
