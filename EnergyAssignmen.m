%% Energy Assignment
clear
clc
R1 = 0.112;
R2d = 0.317;
X1 = j*1.364;
X2d = j*1.32
Xm = j*45.8
% Rc = 0;
N = 94/100; %efficiency 
phi = 0;
Pw = 500e3; %kw
fs = 50;
P = 4;
s =  3.35/100;
parr =@(r1, r2) r1*r2/(r1+r2)
pot = 500e3;

%%
Pin =@(V, I, phi) sqrt(3)*V*I*cos(phi)
Pout =@(T, wm) T*wm
Pdev = @(I, R, s) 3*i^2*R*(1-s)/s
Pag = @(I, R, s) 3*I^2*R/s;

V1 = 2400/sqrt(3)
%%
Pinc = 500e3*0.94;

%% calculation of I1, input current
I1 = V1/(R1+X1+parr(Xm, X2d+R2d/s))
I1mag = abs(I1)
I1phase = radtodeg(angle(I1))
PinC = Pinc*I1';
%% calculation of I2' the load current.
I2d = I1 * Xm /(R2d/s + Xm+X2d)
abs(I2d)
radtodeg(angle(I2d))

%%
pag = Pag(I2d, R2d, s)
abs(pag)
radtodeg(angle(pag))

%%
pdev = Pdev(I2d, R2d, s)
abs(pag)
radtodeg(angle(pag))

ml = (1-N)*Pw - pag - pdev
abs(ml)
radtodeg(angle(ml))

cl = Pw/0.94 - Pinc - ml %copper losses.
abs(cl)
radtodeg(angle(cl))

cl/3

%%
Vp = 2400/sqrt(3)
Ppin =( 500e3 /0.94 )*1/3;
%after finding Rc 
%refind I1
%Pin = V*I'
%pin must be 0.94 * 500e3;
syms rc
z1 = R1 + X1;
z2 = parr(rc, Xm);
z3 = R2d+X2d+R2d/s;
Rt = z1 + parr(z1, z2)
Rc = solve(Rt==Vp^2 /Ppin, rc)
eval(Rc)

%% B
s = -3.2/100;
RC = 400;
Z1 = R1 + X1;
Z2 = parr(RC, Xm);
Z3 = X2+R2d+R2d/s;

