function dx=model(t,x,p)
global Tsp2 Tsp1 Tspan thetad lamd ald cfnom Tfbias tol Qcbias
Q=2.5;
Tfbias=29;
cfnom=0.6;
Qcbias=2;
Tcf=25;
Vc=100;
E=10.1;
K0=2e3;
R=1.98e-3;
ro=850;
Cp=1.35e-4;
deltaHr=-35;
roc=1000;
Cpc=1e-3;
A=900;
U=4e-5;
thetad=0.1;
lamd=1/thetad;
ald=0.05;


yc=p(1);
yi=p(2);
V1=p(3);
V2=p(4);
V=[V1 V2];
Kc1=p(5);
Ki1=p(6);
Kc2=p(7);
Ki2=p(8);


c1=x(1);
T1=x(2);
c2=x(3);
T2=x(4);
Tc1=x(5);
Tc2=x(6);
I1= x(7);
I2= x(8);
J= x(9);
Tc1limit=x(10);
Tc2limit=x(11);
Qclimit=x(12);
TFUlimit=x(13);
TFLlimit=x(14);
Qcu=x(15);

%P1=yo*(Tsp1-T1)+(1-yo)*(Tsp2-T2);
%P2=(1-yo)*(Tsp1-T1)+(yo)*(Tsp2-T2);
P1=(Tsp1-T1);
P2=(Tsp2-T2);
Qc=Qcbias-(1-yi)*(Kc1*P1+Ki1*I1)-yi*(Kc2*P2+Ki2*I2);
Tf=Tfbias+yi*(Kc1*P1+Ki1*I1)+(1-yi)*(Kc2*P2+Ki2*I2);
%
cf=cfnom+ald*(exp(-lamd*t)-1);

Tinc1=Tcf*yc+(1-yc)*Tc2;
Tinc2=Tc1*yc+(1-yc)*Tcf;


theta=V/Q;
beta=deltaHr/(ro*Cp);
alpha=U.*A./(ro.*Cp.*V);
alphac=U*A/(roc*Cpc*Vc);
ra1=-K0*c1*exp(-E/(R*(T1+273.15)));
ra2=-K0*c2*exp(-E/(R*(T2+273.15)));




% eqution
%              cdot=  dc/dt
%              Tdot=  dT/dt

%first reactor
c1dot=(cf-c1)/(theta(1))+ra1;
T1dot=(Tf-T1)/(theta(1))+beta*ra1-alpha(1)*(T1-Tc1);

%second reactor
c2dot=(c1-c2)/(theta(2))+ra2;
T2dot=(T1-T2)/(theta(2))+beta*ra2-alpha(2)*(T2-Tc2);
%for jaket 
Tc1dot=Qc*(Tinc1-Tc1)/(Vc)+alphac*(T1-Tc1);
Tc2dot=Qc*(Tinc2-Tc2)/(Vc)+alphac*(T2-Tc2);




 dI1=P1;
 dI2=P2;
 dJ=1/Tspan*((Tsp2-T2)^2+(Tsp1-T1)^2);
 dTc1f=max(0,-(Tc1-Tcf));
 dTc2f=max(0,-(Tc2-Tcf));
 dQcl=max(0,-Qc);
 dTFUl=max(0,Tf-60);
 dTFLl=max(0,-(Tf));
 dQcu=max(0,Qc-8);
dx=[c1dot;T1dot;c2dot;T2dot;Tc1dot;Tc2dot;dI1;dI2;dJ;dTc1f;dTc2f;dQcl;dTFUl;dTFLl;dQcu];
end