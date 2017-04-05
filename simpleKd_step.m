function dc = simpleKd_step(t,c,constn)

% This function uses the ode45 function to intergrate chemical equations
% for cell volume change simulations. It is used by run_simpleKd.m

A=c(1);
B=c(2);
C=c(3);
k_on=constn(1);
k_off=constn(2);


stoiA=constn(3);
stoiB=constn(4);

dc = [stoiA*k_off*C  - stoiA*k_on*(A^stoiA)*(B^stoiB);
      stoiB*k_off*C  - stoiB*k_on*(A^stoiA)*(B^stoiB);
      - k_off*C      + k_on*(A^stoiA)*(B^stoiB)];
