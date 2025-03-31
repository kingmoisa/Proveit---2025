clear all; close all; clc;
s = tf("s");
G1s = 1 / ((s + 0.01) * (s + 0.1));
% figure; bode(G1s); grid on;

u0=20; g=9.81; m=1000; f=0.015; Theta=0;
rho=1.202; A=1; Cd=0.5; uw=2;
Tau=(m/(rho*A*Cd*(u0+uw)));
K = 1 / (rho * A * Cd * (u0 + uw));

G2s = K / (Tau * s + 1);

Gps = G1s * G2s;
Gps = minreal(Gps);
[numGps, denGps] = tfdata(Gps, 'v');
% step(Gps);grid on;

tt = 5;
suprareglare = 0.02;
% amortizare = abs(log(suprareglare)) / (sqrt(pi^2 + (log(suprareglare))^2));
amortizare = 0.7;
wn = 4 / (amortizare * tt);
i = sqrt(-1);

% p1 = -amortizare * wn + i * wn * sqrt(1-amortizare^2);
% p2 = -amortizare * wn - i * wn * sqrt(1-amortizare^2);
% p3 = -10 * amortizare * wn;
% [A,B,C,D] = tf2ss(numGps, denGps);
% k = place(A,B,[p1  p2  p3]);
% Acl = A - B * k;
% kr=-1*(1/(C*inv(A-B*k)*B));
% sss = ss(Acl, B*kr, C, D);
% step(70*sss);
p1 = -amortizare * wn + i * wn * sqrt(1-amortizare^2);
p2 = -amortizare * wn - i * wn * sqrt(1-amortizare^2);
p3 = -5 * wn;
[A,B,C,D] = tf2ss(numGps, denGps);
k = place(A,B,[p1  p2  p3]);
Acl = A - B * k;
kr = -1*(1/(C*inv(A-B*k)*B));
sss = ss(Acl, B*kr, C, D)

% G0s = wn^2 / (s^2 + 2 * amortizare * wn * s + wn ^ 2);

[numG0s, denG0s] = tfdata(G0s, 'v');
Grs = (1 / Gps) *  (G0s/ (1 - G0s));
Grs = minreal(Grs);
[numGrs, denGrs] = tfdata(Grs, 'v');

Ku = 0.25782;
Tu = 134.449;
C=10^-1;
D = 10^0;
Kp = Ku * 0.60;
Ki = 1.2 * Ku / Tu * C;
Kd = 0.075 * Ku / Tu * D;