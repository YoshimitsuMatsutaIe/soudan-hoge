
#include "A.h"

#include "A_0_0.h"
#include "A_0_1.h"
#include "A_0_2.h"
#include "A_0_3.h"
#include "A_0_4.h"
#include "A_0_5.h"
#include "A_0_6.h"
#include "A_0_7.h"
#include "A_0_8.h"

#include "A_1_0.h"
#include "A_1_1.h"
#include "A_1_2.h"
#include "A_1_3.h"
#include "A_1_4.h"
#include "A_1_5.h"
#include "A_1_6.h"
#include "A_1_7.h"
#include "A_1_8.h"

#include "A_2_0.h"
#include "A_2_1.h"
#include "A_2_2.h"
#include "A_2_3.h"
#include "A_2_4.h"
#include "A_2_5.h"
#include "A_2_6.h"
#include "A_2_7.h"
#include "A_2_8.h"

#include "A_3_0.h"
#include "A_3_1.h"
#include "A_3_2.h"
#include "A_3_3.h"
#include "A_3_4.h"
#include "A_3_5.h"
#include "A_3_6.h"
#include "A_3_7.h"
#include "A_3_8.h"

#include "A_4_0.h"
#include "A_4_1.h"
#include "A_4_2.h"
#include "A_4_3.h"
#include "A_4_4.h"
#include "A_4_5.h"
#include "A_4_6.h"
#include "A_4_7.h"
#include "A_4_8.h"

#include "A_5_0.h"
#include "A_5_1.h"
#include "A_5_2.h"
#include "A_5_3.h"
#include "A_5_4.h"
#include "A_5_5.h"
#include "A_5_6.h"
#include "A_5_7.h"
#include "A_5_8.h"

#include "A_6_0.h"
#include "A_6_1.h"
#include "A_6_2.h"
#include "A_6_3.h"
#include "A_6_4.h"
#include "A_6_5.h"
#include "A_6_6.h"
#include "A_6_7.h"
#include "A_6_8.h"

#include "A_7_0.h"
#include "A_7_1.h"
#include "A_7_2.h"
#include "A_7_3.h"
#include "A_7_4.h"
#include "A_7_5.h"
#include "A_7_6.h"
#include "A_7_7.h"
#include "A_7_8.h"

#include "A_8_0.h"
#include "A_8_1.h"
#include "A_8_2.h"
#include "A_8_3.h"
#include "A_8_4.h"
#include "A_8_5.h"
#include "A_8_6.h"
#include "A_8_7.h"
#include "A_8_8.h"

void A(double h1, double h2, double h3, double l1, double l1_dot, double l2, double l2_dot, double l3, double l3_dot, double *out){
    out[0] = A_0_0();
    out[1] = A_1_0();
    out[2] = A_2_0();
    out[3] = A_3_0(h1, h2, h3, l1, l1_dot, l2, l2_dot, l3, l3_dot);
    out[4] = A_4_0(h1, h2, h3, l1, l1_dot, l2, l2_dot, l3, l3_dot);
    out[5] = A_5_0(h1, h2, h3, l1, l1_dot, l2, l2_dot, l3, l3_dot);
    out[6] = A_6_0();
    out[7] = A_7_0();
    out[8] = A_8_0();
    out[9] = A_0_1();
    out[10] = A_1_1();
    out[11] = A_2_1();
    out[12] = A_3_1(h1, h2, h3, l1, l1_dot, l2, l2_dot, l3, l3_dot);
    out[13] = A_4_1(h1, h2, h3, l1, l1_dot, l2, l2_dot, l3, l3_dot);
    out[14] = A_5_1(h1, h2, h3, l1, l1_dot, l2, l2_dot, l3, l3_dot);
    out[15] = A_6_1();
    out[16] = A_7_1();
    out[17] = A_8_1();
    out[18] = A_0_2();
    out[19] = A_1_2();
    out[20] = A_2_2();
    out[21] = A_3_2(h1, h2, h3, l1, l1_dot, l2, l2_dot, l3, l3_dot);
    out[22] = A_4_2(h1, h2, h3, l1, l1_dot, l2, l2_dot, l3, l3_dot);
    out[23] = A_5_2(h1, h2, h3, l1, l1_dot, l2, l2_dot, l3, l3_dot);
    out[24] = A_6_2();
    out[25] = A_7_2();
    out[26] = A_8_2();
    out[27] = A_0_3();
    out[28] = A_1_3();
    out[29] = A_2_3();
    out[30] = A_3_3(l1, l1_dot, l2, l2_dot, l3, l3_dot);
    out[31] = A_4_3(l1, l1_dot, l2, l2_dot, l3, l3_dot);
    out[32] = A_5_3(l1, l1_dot, l2, l2_dot, l3, l3_dot);
    out[33] = A_6_3();
    out[34] = A_7_3();
    out[35] = A_8_3();
    out[36] = A_0_4();
    out[37] = A_1_4();
    out[38] = A_2_4();
    out[39] = A_3_4(l1, l1_dot, l2, l2_dot, l3, l3_dot);
    out[40] = A_4_4(l1, l1_dot, l2, l2_dot, l3, l3_dot);
    out[41] = A_5_4(l1, l1_dot, l2, l2_dot, l3, l3_dot);
    out[42] = A_6_4();
    out[43] = A_7_4();
    out[44] = A_8_4();
    out[45] = A_0_5();
    out[46] = A_1_5();
    out[47] = A_2_5();
    out[48] = A_3_5(l1, l1_dot, l2, l2_dot, l3, l3_dot);
    out[49] = A_4_5(l1, l1_dot, l2, l2_dot, l3, l3_dot);
    out[50] = A_5_5(l1, l1_dot, l2, l2_dot, l3, l3_dot);
    out[51] = A_6_5();
    out[52] = A_7_5();
    out[53] = A_8_5();
    out[54] = A_0_6();
    out[55] = A_1_6();
    out[56] = A_2_6();
    out[57] = A_3_6(l1, l2, l3);
    out[58] = A_4_6(l1, l2, l3);
    out[59] = A_5_6(l1, l2, l3);
    out[60] = A_6_6();
    out[61] = A_7_6();
    out[62] = A_8_6();
    out[63] = A_0_7();
    out[64] = A_1_7();
    out[65] = A_2_7();
    out[66] = A_3_7(l1, l2, l3);
    out[67] = A_4_7(l1, l2, l3);
    out[68] = A_5_7(l1, l2, l3);
    out[69] = A_6_7();
    out[70] = A_7_7();
    out[71] = A_8_7();
    out[72] = A_0_8();
    out[73] = A_1_8();
    out[74] = A_2_8();
    out[75] = A_3_8(l1, l2, l3);
    out[76] = A_4_8(l1, l2, l3);
    out[77] = A_5_8(l1, l2, l3);
    out[78] = A_6_8();
    out[79] = A_7_8();
    out[80] = A_8_8();
}
