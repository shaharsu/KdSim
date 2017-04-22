function sse=sseRun_simpleKd(x0,start,exp_V,exp_cR_avg,exp_cG_avg)
% Usage:
% 
% [sse,ks_on,ks_off,R_i,G_i,sim_cR,sim_cG,sim_cF,dV,A_i,B_i,C_i]=sseRun_simpleKd(k_on, k_off)
% 
% This function calculates the sum of square errors between a series of
% n_run simulations and experimental chi data.
% It is to be used by the minimizeRun_simpleKd.m function.
% currently contains all data used in the original paper.


n_runs  = 20; % number of simulations

A=start(1);
B=start(2);
C=start(3);
dV=start(4);
k_on=10^(x0);
k_off=start(6);
E_C = start(7); % FRET efficiency of complex C
stoiA=start(8);
stoiB=start(9);

% assign vectors for simulation results
sim_cR = zeros(n_runs,6);
sim_cG = zeros(n_runs,6);
sim_cF = zeros(n_runs,6);

%Generate result sets


for i=1:size(exp_V,2)
    for j=1:n_runs
        [A_in,B_in]=deal(-1);
        while (A_in<0)
            A_in = (A*randn(1)+5)*1e-6; % A from random variable
        end
        while (B_in/A_in<0.5||B_in/A_in>2)
            B_in = (B*randn(1)+6)*1e-6; % B from random variables
        end
        dV_i = exp_V(i)+randn()*0.01; % random dV
        start=[A_in B_in 0 dV_i k_on k_off E_C, stoiA, stoiB]; % random dV
        %[tmp1(i),tmp2(i)]=run_simpleKd(start)
        [~,~,sim_cR(j,i), sim_cG(j,i),~,~,~,~,~,~,~] = run_simpleKd(start);
    end
end

sim_cR_avg = mean(sim_cR,1,'omitnan');
sim_cG_avg = mean(sim_cG,1,'omitnan');
sse = sum((sim_cG_avg-exp_cG_avg).^2) + sum((sim_cR_avg-exp_cR_avg).^2);


end

