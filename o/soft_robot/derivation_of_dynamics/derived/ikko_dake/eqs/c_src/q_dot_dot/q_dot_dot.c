#include "q_dot_dot.h"



#include "q_dot_dot_0.h"
#include "q_dot_dot_1.h"
#include "q_dot_dot_2.h"


void q_dot_dot(double h1, double h2, double h3, double l1, double l1_dot, double l2, double l2_dot, double l3, double l3_dot, double tau1, double tau2, double tau3, double *out) {
    out[0] = q_dot_dot_0(h1, h2, h3, l1, l1_dot, l2, l2_dot, l3, l3_dot, tau1, tau2, tau3);
    out[1] = q_dot_dot_1(h1, h2, h3, l1, l1_dot, l2, l2_dot, l3, l3_dot, tau1, tau2, tau3);
    out[2] = q_dot_dot_2(h1, h2, h3, l1, l1_dot, l2, l2_dot, l3, l3_dot, tau1, tau2, tau3);
}
