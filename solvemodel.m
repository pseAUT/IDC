function [t,x]=solvemodel(p)
global  Tspan
global x0 V yc
V=[p(3),p(4)];
yc=p(1);
% % c10=0.0675517011982044;
% % T10=173.428991859559;
% % c20=0.00760744461288107;
% % T20=173.417075427106;
% % Tc10=30.2023500416110;
% % Tc20=27.6242704888879;
% % I10= 0;
% % I20= 0;
% % J0= 0;
% % Tc1limit0=0;
% % Tc2limit0=0;
% % Qclimit=0;
% % TFUlimit=0;
% % TFLlimit=0;
% % Qcu=0;

% x0=[c10;
% T10;
% c20;
% T20;
% Tc10;
% Tc20;
% I10;
% I20;
% J0;
% Tc1limit0;
% Tc2limit0;
% Qclimit;
% TFUlimit;
% TFLlimit;
% Qcu];
x0= fsolve(@INIG,800*ones(1,6));
[t,x]=ode23s(@(t,x)model(t,x,p),[0 Tspan],[x0';0;0;0;0;0;0;0;0;0]);





end