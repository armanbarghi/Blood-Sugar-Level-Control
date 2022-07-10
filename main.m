clear; clc;
Gb = 4.5;
Ib = 15;
P1 = 0;
P2 = 0.025;
P3 = 0.000013;
V = 12;
n = 5.54/60;
syms G X I U D
F = [ -P1*G-X*(G+Gb)+D;
      -P2*X+P3*I;
      -n*(I+Ib)+U/V ];
y = G;

St = [G;X;I];
In = [U;D];

%% Equilibrium State
% sol_s = solve(subs(F,{D},{0})==[0;0;0], [G,X,I]);
% Ge = sol_s.G; Ie = sol_s.I; Xe = sol_s.X;
St_e = [0;0;0];
In_e = [n*Ib*V;0];
Fe = subs(subs(F,In,In_e),St,St_e);

%% Linearize
% f(x) = f(xe,ue) + (grad(f,x)|x=xe,u=ue)*(x-xe) + (grad(f,u)|x=xe,u=ue)*(u-ue)
grad_F_St = [ diff(F(1),St(1)), diff(F(1),St(2)), diff(F(1),St(3));
              diff(F(2),St(1)), diff(F(2),St(2)), diff(F(2),St(3));
              diff(F(3),St(1)), diff(F(3),St(2)), diff(F(3),St(3)) ];
grad_F_In = [ diff(F,In(1)) , diff(F,In(2)) ];
A = subs(subs(grad_F_St,In,In_e),St,St_e)
B = subs(subs(grad_F_In,In,In_e),St,St_e)
C = [1,0,0]
D = 0

%% Transfer Function
syms s
n = size(A,1);
G = C*inv(s*eye(n)-A)*B;
ExpFun = matlabFunction(simplifyFraction(G,'Expand',true));
ExpFun = str2func(regexprep(func2str(ExpFun), '\.([/^\\*])', '$1'));
G = tf(ExpFun(tf('s')));
for i = 1:length(G)
    [num,den] = tfdata(G(i));
    G(i) = tf(num{1}/den{1}(1),den{1}/den{1}(1));
end
G