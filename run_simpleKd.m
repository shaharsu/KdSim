function [R_i,G_i,dR,dG,dF,dV,A_i,B_i,C_i,k_on,k_off,t_tot,A,B,C,R,G,F] = run_simpleKd(start)
%
% run_simpleKd is the main function that generates cell volume change data.
% It can be used on its own, but is also used by KdSimGUI.m and
% minimizeRun_simpleKd.m to generates chi data.
%
% The optional input vairable "start" can supply initial conditions:
% start = [A B C dV k_on k_off E_C stoiA stoiB]
% you can also just supply k_on and k_off, or nothing at which uses default
% settings shown below.

if ~exist('start','var')
    A(1)=5e-6;          % donor
    B(1)=4.5e-6;      %*(abs(rand())+0.5);    % acceptor
    C(1)=0;             % FRET complexS
    dV=0.75; %rand()+0.5;
    k_on=1e5;
    k_off=1e-1;
    E_C = 0.3;          % FRET efficiency of complex C
    showfig=1;
    stoiA=1;
    stoiB=1;
elseif max(size(start))==2
    A(1)=3e-6;          % donor
    B(1)=A(1)*0.9;      %*(abs(rand())+0.5);    % acceptor
    C(1)=0;             % FRET complexS
    dV=0.75; %rand()+0.5;
    k_on=start(1);
    k_off=start(2);
    E_C = 0.3; % FRET efficiency of complex C
    stoiA=4;
    stoiB=1;
    showfig=1;
else
    A=start(1);
    B=start(2);
    C=start(3);
    dV=start(4);
    k_on=start(5);
    k_off=start(6);
    E_C = start(7); % FRET efficiency of complex C
    stoiA=start(8);
    stoiB=start(9);
    showfig=0;
end

A_i=A;
B_i=B;
C_i=C;


% Check if timescale is too long for rates (expedites intergration)

if log10(k_on)/(stoiA+stoiB-1)<9
    tspan=0:0.1:100;
else
    tspan=0:1e-5:0.01;
end

% Eq. constants
c=[A,B,C];
K_D=k_off/k_on;
params = [k_on, k_off, stoiA, stoiB];
opts_1 = odeset('RelTol',1e-8,'AbsTol',1e-9);%,'outputfcn',@odeplot,'events',@checkConv);

%--------------------------------
% integrate the function simpleKd_step using ode45
%--------------------------------

[t,c]=ode45(@(t,c)simpleKd_step(t,c,params),tspan,c,opts_1);
%volume change
c2=c(end,:)/dV;
[t2,c2]=ode45(@(t,c)simpleKd_step(t,c,params),tspan,c2,opts_1);

%total signal
t_tot=[t;t(end)+t2];
l=[(t+1)./(t+1);(t2+1)-(t2+1)+dV];
c_tot=[c;c2];
A=c_tot(:,1);
B=c_tot(:,2);
C=c_tot(:,3);

%-------------------------------
% calculate fluorescence from concentrations
%-------------------------------

% fluorescence parameters
S_g = 80e3;   % cross section of acceptor at 475 nm
S_r = 10e3;   % cross section of donor at 475 nm
QY_A = 0.8;   % quantum yield donor
QY_B = 0.3;   % quantum yield acceptor

R = l.*(B.*S_r.*QY_B + ... % Red cross-excitation from acceptor under blue light
    S_r.*QY_B.*C.*(1-E_C) + ...  % Red cross-excitation from acceptor in complex under blue light
    S_g.*QY_B.*C.*(E_C));% + ...  % Red FRET from complex


G = l.*((1-E_C).*C.*S_g.*QY_A + ... % Green from complex
    A.*S_g.*QY_A); % + ... % Green from donor


i0=size(t,1)+2;
F = R./(R+G);
R_i=R(i0-2);
G_i=G(i0-2);
F_i=F(i0-2);
dF=(F(end)-F_i)/F_i;

if R_i==0
    dR=0;
else
    dR=(R(end)-R_i)/R_i;
end

if G_i==0
    dG=0;
else
    dG=(G(end)-G_i)/G_i;
end

%-------------------------------
% Draw some pictures, if showfig==1
%-------------------------------

if showfig==1
    figure('position',[0 100 800 500]); % [left bottom width height]
    t_ini=1;
    t_fin=max(size(t_tot));
    
    sh(1)=subplot(2,2,1);
    plot(t_tot,[A(t_ini:t_fin),B(t_ini:t_fin),C(t_ini:t_fin)],'linewidth',0.5,'marker','.')
    line([t_tot(i0-1-t_ini) t_tot(i0-1-t_ini)],get(gca,'ylim'),'color','black','linestyle','--')
    legend('Donor','Acceptor','Complex');
    xlabel('t (s)');
    ylabel('C (M)');
    
    sh(2)=subplot(2,2,2);
    hax1=plot(t_tot,G(t_ini:t_fin),'marker','.','color','g','linewidth',0.5);
    line([t_tot(i0-1-t_ini) t_tot(i0-1-t_ini)],get(gca,'ylim'),'color','black','linestyle','--')
    xlabel('t (s)');
    ylabel('F_{green}');
    
    sh(3)=subplot(2,2,4);
    hax2=plot(t_tot,R(t_ini:t_fin),'marker','.','color','r','linewidth',0.5);
    line([t_tot(i0-1-t_ini) t_tot(i0-1-t_ini)],get(gca,'ylim'),'color','black','linestyle','--')
    xlabel('t (s)');
    ylabel('F_{red}');
    
    sh(4)=subplot(2,2,3);
    plot(t_tot,F(t_ini:t_fin)./F(t_ini),'linewidth',0.5,'marker','.');
    line([t_tot(i0-1-t_ini) t_tot(i0-1-t_ini)],get(gca,'ylim'),'color','black','linestyle','--')
    xlabel('t (s)');
    ylabel('normalized E_{FRET}');
    line1=sprintf('K_D = %1.e\nk_{on} = %.1e\nd_{off} = %.1e\ndR = %.1e\ndG = %1.e\ndV = %.2f\n',[K_D,k_on,k_off,dR,dG,dV]);
    textpos=[get(gca,'Xlim') get(gca,'Ylim')];
    text(textpos(1)+textpos(2)*0.05,textpos(4)-textpos(3)*0.1,line1);
    linkaxes(sh,'x');
    
end


end
