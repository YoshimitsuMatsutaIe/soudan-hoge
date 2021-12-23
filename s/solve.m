
syms LAMBDA CONTACT_E TOTAL CONTACT_I CONTACT_S1 ALPHA MU RATE_W GAMMA...
    THETA1 THETA2 SPAN1 SPAN2 SPAN12 S E I R S1

eq1 = LAMBDA - (CONTACT_E*S*E/TOTAL+ CONTACT_I*S*I/TOTAL) - S*(SPAN1*THETA1+MU);
eq2 =  (CONTACT_E*S*E/TOTAL+ CONTACT_I*S*I/TOTAL+CONTACT_S1*S1*I/TOTAL+CONTACT_S1*S1*E/TOTAL) - E*(MU+ALPHA);
eq3 =  ALPHA*E - I*(RATE_W + MU + GAMMA);
eq4 =  S1*SPAN2*THETA2*SPAN12 + GAMMA*I - R*MU;
eq5 =  S*SPAN1*THETA1- S1*(SPAN2*THETA2*SPAN12 + MU)-(CONTACT_S1*S1*I/TOTAL)-(CONTACT_S1*S1*E/TOTAL);

sol = solve(eq1==0, eq2==0, eq3==0, eq4==0, eq5==0, S, E, I, R, S1);


matlabFunction(sol.S,'file','endemic_sol_S.m')
matlabFunction(sol.E,'file','endemic_sol_E.m')
matlabFunction(sol.I,'file','endemic_sol_I.m')
matlabFunction(sol.R,'file','endemic_sol_R.m')
matlabFunction(sol.S1,'file','endemic_sol_S1.m')