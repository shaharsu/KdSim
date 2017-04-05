function sse=sseRun_simpleKd(x0,start)
% Usage:
% 
% [sse,ks_on,ks_off,R_i,G_i,sim_cR,sim_cG,sim_cF,dV,A_i,B_i,C_i]=sseRun_simpleKd(k_on, k_off)
% 
% This function calculates the sum of square errors between a series of
% n_run simulations and experimental chi data.
% It is to be used by the minimizeRun_simpleKd.m function.
% currently contains all data used in the original paper.

n_runs  = 10; % number of simulations

A=start(1);
B=start(2);
C=start(3);
dV=start(4);
k_on=10^(x0);
k_off=start(6);
E_C = start(7); % FRET efficiency of complex C
stoiA=start(8);
stoiB=start(9);

% simulated results, n_runs=10, A=B=5e-6, E_C=0.3, stoiA=stoiB=1, kD=1e-6
exp_V = [1.25218246804198,1.06680002836414,1.00416259201074,0.921035706316025,0.831432066656626,0.740767351501897];
exp_cR_avg = [-0.0259797841317013,-0.00764197576232570,-0.000453467419212500,0.00915410694516710,0.0195633606342439,0.0320401485999640];
exp_cG_avg = [0.0119309536884247,0.00313972560103450,0.000208745204338100,-0.00359946227257470,-0.00845718314999690,-0.0128045195638633];

% simulated results, n_runs=10, A=B=5e-6, E_C=0.3, stoiA=3, stoiB=1, kD=1e-13
% exp_V=[1.25100004660845,1.06617985921542,0.998893070682346,0.924901831203719,0.828939779971891,0.743150021714604];
% exp_cR_avg=[-0.00491184380761650,-0.00164481024675500,7.50421115837000e-05,0.00107437377672380,0.00771052782158320,0.0234033803577589];
% exp_cG_avg=[0.00417805041518350,0.00145504674972940,-5.16996065732000e-05,-0.000797673850917900,-0.00896914782424070,-0.0118327539534502];

% simulated results, n_runs=10, A=B=5e-6, E_C=0.3, stoiA=2, stoiB=2, kD=1e-17
% exp_V=[1.24668685758402,1.07110002256849,1.00040929215816,0.921277060211530,0.828879124781285,0.740052788997494];
% exp_cR_avg=[-0.0168423016399495,-0.00343790434504540,-5.14426207224000e-05,0.00677686107270260,0.0114220888158051,0.0152602020335052];
% exp_cG_avg=[0.0451210874800126,0.0115617323141159,-0.000533821091849900,-0.0109816721430914,-0.0357135132372920,-0.0546274173713987];

% Experimental results - AcGFP1:mCherry
% exp_V   = [1.25 1.07 1 0.92 0.83 0.74]; % Experimental average cell volume changes
% exp_cR_avg = [-0.084 -0.03 0 0.04 0.1 0.16];
% exp_cG_avg = [0.16 0.02 0 -0.03 -0.04 -0.1];

% Experimental results - GAPDH-mEGFP:PGK-mCherry
% exp_V   = [1.25 1.07 1 0.92 0.83 0.74]; % Experimental average cell volume changes
% exp_cR_avg = [-0.03832983452888 -0.02374839764371 0 0.0191802499543 0.1257409433196 0.1538579687254];
% exp_cG_avg = [0.0337973202795 -0.00127996691427 0 -0.01719731500701 -0.02467624195907 -0.06346443955339];

% assign vectors for simulation results
sim_cR = zeros(n_runs,6);
sim_cG = zeros(n_runs,6);
sim_cF = zeros(n_runs,6);

%Generate result sets

for i=1:6
    for j=1:n_runs
        [A_in,B_in]=deal(-1);
        while (A_in<0)
            A_in = (A*randn(1)+6)*1e-6; % A from random variable
        end
        while (B_in<0)
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

