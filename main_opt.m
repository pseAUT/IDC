
clear all
close all
clc
global Tsp2 Tsp1 Tspan thetad lamd ald cfnom Tfbias tol Qcbias
global x0

format long

tol=1e-10;
Tspan=3000;
Tsp1=0.9*173;
Tsp2=0.9*173;

nvars=8;
A=[];
b=[];
Aeq=[];
beq=[];
LB=[0 ;0 ;720 ;720 ;0.00 ;0.0 ;0 ;0.0];
UB=[1 ;1 ;1080 ;1080 ;50 ;50 ;50 ;50];
nonlcon=@cons;
Intcon=[1;2];
options=[];

tic
[p,fval,flag,output,population,scores] = ga(@obj,nvars,A,b,Aeq,beq,LB,UB,nonlcon,Intcon,options);
% [p, OBJ, INFO]=ga(@obj,2,[],[],[],[],LB,UB,@cons)
toc;
flag

yc=p(1);
yi=p(2);
%yo=p(3);
V1=p(3);
V2=p(4);
V=[V1 V2];
Kc1=p(5);
Ki1=p(6);
Kc2=p(7);
Ki2=p(8);

%Plot results to verify performance & stability
[t,x]=solvemodel(p);
cf=cfnom+ald*(exp(-lamd*t)-1);


c1=x(:,1);
T1=x(:,2);
c2=x(:,3);
T2=x(:,4);
Tc1=x(:,5);
Tc2=x(:,6);
I1= x(:,7);
I2= x(:,8);
J= x(:,9);

P1=(Tsp1-T1);
P2=(Tsp2-T2);
Tf=Tfbias+yi*(Kc1*P1+Ki1*I1)+(1-yi)*(Kc2*P2+Ki2*I2);
Qc=Qcbias-(1-yi)*(Kc1*P1+Ki1*I1)-yi*(Kc2*P2+Ki2*I2);



  subplot(3,4,1)
    
    plot(t,cf,'r','LineWidth',2)
    ylabel('cf')
   xlabel('time')
    grid on
    
    subplot(3,4,2)
    plot(t,c1,'LineWidth',2)
    ylabel('c1')
    xlabel('time')
    grid on
    
    subplot(3,4,3)
    plot(t,T1,'LineWidth',2)
    ylabel('T1')
    xlabel('time')
    grid on
    
    subplot(3,4,4)
    plot(t,c2,'LineWidth',2)
    ylabel('c2')
    xlabel('time')
    grid on
    
      subplot(3,4,5)
    plot(t,T2,'LineWidth',2)
    ylabel('T2')
    xlabel('time')
    grid on
    
      subplot(3,4,6)
    plot(t,Tc1,'LineWidth',2)
    ylabel('Tc1')
    xlabel('time')
    grid on
    
      subplot(3,4,7)
    plot(t,Tc2,'LineWidth',2)
    ylabel('Tc2')
    xlabel('time')
    grid on
      subplot(3,4,8)
    plot(t,Tf,'LineWidth',2)
    ylabel('Tf')
    xlabel('time')
    grid on
    subplot(3,4,9)
    plot(t,Qc,'LineWidth',2)
    ylabel('QC')
    xlabel('time')
    grid on