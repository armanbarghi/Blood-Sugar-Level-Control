clear; clc;
Gb = 4.5;
Ib = 15;
P1 = 0;
P2 = 0.025;
P3 = 0.000013;
V = 12;
n = 5.54/60;
syms G_dot X_dot I_dot G X I D U
F = [ -P1*G-X*(G+Gb)+D;
      -P2*X+P3*I;
      -n*(I+Ib)+U/V ];
y = G;

%% Equilibrium State
% states = solve(subs(F,{D},{0})==[0;0;0], [G,X,I]);
% Ge = states.G; Ie = states.I; Xe = states.X;
Ge = Gb; Xe = 0; Ie = 0; Ue = n*Ib*V;
Fe = subs(F,{G,X,I,U},{Ge,Xe,Ie,Ue});

%% Linearize
% f(x) = f(x0) + transpose(grad(f,x)|x=x0)*(x-x0)
grad_F = [ diff(F(1),G), diff(F(1),X), diff(F(1),I);
           diff(F(2),G), diff(F(2),X), diff(F(2),I);
           diff(F(3),G), diff(F(3),X), diff(F(3),I) ];
lin_F = Fe + subs(grad_F,{G,X,I,U},{Ge,Xe,Ie,Ue})*([G;X;I]-[Ge;Xe;Ie]);
A = [ diff(lin_F(1),G), diff(lin_F(1),X), diff(lin_F(1),I);
      diff(lin_F(2),G), diff(lin_F(2),X), diff(lin_F(2),I);
      diff(lin_F(3),G), diff(lin_F(3),X), diff(lin_F(3),I) ];
B = [ diff(lin_F(1),U) ; diff(lin_F(2),U) ; diff(lin_F(3),U) ];
vpa(lin_F,3)
vpa(A,3)
vpa(B,3)
D = 0;
C = [1,0,0];

%% Transfer Function
syms s
n = size(A,1);
G = C*inv(s*eye(n)-A)*B
ExpFun = matlabFunction(simplifyFraction(G,'Expand',true));
ExpFun = str2func(regexprep(func2str(ExpFun), '\.([/^\\*])', '$1'));
G = tf(ExpFun(tf('s')));
for i = 1:length(G)
    [num,den] = tfdata(G(i));
    G(i) = tf(num{1}/den{1}(1),den{1}/den{1}(1));
end
G