function [m, R, e] = Kalman_Filter(Y, G, W)
% m(t)     = Estimated State
% R        = Estimation Error
% e        = Forecast Error
N = length(Y);
%% Initialization of the model
m = zeros(1,N);
C = zeros(1,N);

F = 1;
V = 0;

a = zeros(1,N);
R = zeros(1,N);
f = zeros(1,N);
Q = zeros(1,N);
e = zeros(1,N);
%% Iteration 0
m0 = Y(1);
C0 = 1;
%% Iteration 1
a(1) = G*m0;
R(1) = G*C0*G+W;
f(1) = F*a(1);
Q(1) = F*R(1)*F+V;
e(1) = Y(1)-f(1);
m(1) = a(1) + R(1)*F*(1/Q(1))*e(1);
C(1) = R(1) - R(1)*F*(1/Q(1))*F*R(1);
%%
for i = 2:N
    a(i) = G*m(i-1);
    R(i) = G*C(i-1)*G+W;
    f(i) = F*a(i);
    Q(i) = F*R(i)*F+V;
    e(i) = Y(i)-f(i);
    m(i) = a(i) + R(i)*F*(1/Q(i))*e(i);
    C(i) = R(i) - R(i)*F*(1/Q(i))*F*R(i);
end
R = m-Y;
end