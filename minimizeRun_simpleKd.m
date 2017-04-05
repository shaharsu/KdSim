function [landscape,x,fval,eflag,output,manymins]=minimizeRun_simpleKd(start)

% This function performs a global minimization for a set stoichiometry. The
% function generates starting points using the multistart function from a
% range determined by the stoichiometry. multistart is run in parallel to
% make it more efficient.
%
% details of the input parameter, start, is specified below
% Generally this function is run within a for loop trying out several
% different stoichiometry combinations. Running the following lines will
% perform this analysis:
%
% stoi=[1 1;1 2;1 3;1 4;1 5;2 1;3 1;4 1;5 1;2 2;2 3;3 2]
% A=cell(1,12);
% for i=1:12
%    display(sprintf('!!! on run %i',i));
%    A{i}=minimizeRun_simpleKd([5 5 0 0.75 1e5 1e-1 0.6 stoi(i,1) stoi(i,2)]);
%    B=A{i}
%    dlmwrite(strcat(num2str(stoi(i,1)),'to',num2str(stoi(i,2)),'.dat')
% end 
%
% The output .dat files contain two columns: The log(k_on) and the sse (sum
% of square errors)for each of the points generated by multistart, ordered
% by increasing sse. Remember that
% log(K_D)=log(k_off)-log(k_on)=-(1+log(k_on) if k_off = 0.1
%

if ~exist('start','var')
    start(1)=5;             % [donor] average(uM)
    start(2)=5;             % [acceptor]average(uM
    start(3)=0;             % [complex] average(uM)
    start(4)=0.75;          % deprecated
    start(5)=1e5;           % k_on;
    start(6)=1e-1;          % k_off
    start(7)= 0.6;          % FRET efficiency of complex C
    start(8)=1;             % stoiA
    start(9)=1;             % stoiB
end

x0 = (start(8)+start(9)-1)*5;
lb=(0.5*x0);
ub=(1.2*x0);

%     fitOpts=optimset('display','iter','maxfunevals',4e2,'tolfun',1e-4,'tolX',1e-4,'PlotFcns',@optimplotfval);
%     [x, fval]=fminsearch('sseRunByConc_simpleKd',10,fitOpts);
opts = optimoptions(@fmincon,'Algorithm','interior-point','display','none','steptolerance',1e-4);
problem = createOptimProblem('fmincon','objective',@(x)sseRun_simpleKd(x,start),'x0',x0,'lb',lb,'ub',ub,'options',opts);
ms = MultiStart('UseParallel',true,'Display','iter');
[x,fval,eflag,output,manymins] = run(ms,problem,50);

% globalsearch function may be used by uncommenting this line.
% gs = GlobalSearch('Display','iter','PlotFcn',@gsplotbestf,'FunctionTolerance',1e-3);
%[x,fval, exitflag, output, solutions] = run(gs,problem);

landscape=[manymins.X ;manymins.Fval]';
end
