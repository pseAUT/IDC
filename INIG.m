function F=INIG(x)
global V yc



Q=2.5;
Tf=29;
cf=0.6;
Qc=2;
Tcf=25;
Vc=100;
E=10.1;
K0=2e3;
R=1.98e-3;
ro=850;
%V=900;
Cp=1.35e-4;
deltaHr=-35;
roc=1000;
Cpc=1e-3;
A=900;
U=4e-5;
thetad=0.1;
lamd=1/thetad;
ald=0.05;


c1=x(1);
T1=x(2);
c2=x(3);
T2=x(4);
Tc1=x(5);
Tc2=x(6);

theta=V/Q;
beta=deltaHr/(ro*Cp);
alpha=U.*A./(ro.*Cp.*V);
alphac=U*A/(roc*Cpc*Vc);
ra1=-K0*c1*exp(-E/(R*(T1+273.15)));
ra2=-K0*c2*exp(-E/(R*(T2+273.15)));
Tinc1=Tcf*yc+(1-yc)*Tc2;
Tinc2=Tc1*yc+(1-yc)*Tcf;


%first reactor
c1dot=(cf-c1)/(theta(1))+ra1;
T1dot=(Tf-T1)/(theta(1))+beta*ra1-alpha(1)*(T1-Tc1);

%second reactor
c2dot=(c1-c2)/(theta(2))+ra2;
T2dot=(T1-T2)/(theta(2))+beta*ra2-alpha(2)*(T2-Tc2);
%for jaket 
Tc1dot=Qc*(Tinc1-Tc1)/(Vc)+alphac*(T1-Tc1);
Tc2dot=Qc*(Tinc2-Tc2)/(Vc)+alphac*(T2-Tc2);

F=[c1dot;T1dot;c2dot;T2dot;Tc1dot;Tc2dot];
end