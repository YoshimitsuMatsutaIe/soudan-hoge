#include "invM.h"


#include "invM_0_0.h"
#include "invM_0_1.h"
#include "invM_0_2.h"

#include "invM_1_0.h"
#include "invM_1_1.h"
#include "invM_1_2.h"

#include "invM_2_0.h"
#include "invM_2_1.h"
#include "invM_2_2.h"




void invM(double l1, double l2, double l3, double *out){
    out[0] = invM_0_0(l1, l2, l3);
    out[1] = invM_1_0(l1, l2, l3);
    out[2] = invM_2_0(l1, l2, l3);
    out[3] = invM_0_1(l1, l2, l3);
    out[4] = invM_1_1(l1, l2, l3);
    out[5] = invM_2_1(l1, l2, l3);
    out[6] = invM_0_2(l1, l2, l3);
    out[7] = invM_1_2(l1, l2, l3);
    out[8] = invM_2_2(l1, l2, l3);
}
