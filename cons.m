function [C,Ceq]=cons(p)
global tol


yc=p(1);
yi=p(2);
V1=p(3);
V2=p(4);
V=[V1 V2];
Kc1=p(5);
Ki1=p(6);
Kc2=p(7);
Ki2=p(8);

[t,x]=solvemodel(p);

%Flow lower limit

C(1)=x(end,10)-tol;
C(2)=x(end,11)-tol;
C(3)=x(end,12)-tol;
C(4)=x(end,13)-tol;
C(5)=x(end,14)-tol;
C(6)=x(end,15)-tol;
Ceq=[];
end



