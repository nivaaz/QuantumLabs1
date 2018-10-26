%% Question 1
clear
clc

parr = @(one, two) one*two/(one + two); %parrallel eqn
Vphase = 2400/sqrt(3)
s = 3.35/100;

Pout = 500e3;
R1 = 0.112;
R2 = 0.317;
X1 = j*1.364;
X2 = j*1.32
Xm = j*45.8

Z1 = R1 + X1;
Z2 = Xm;
Z3 = R2/s + X2;

RT = Z1 + parr(Z2, Z3)

I1 = Vphase/RT
radtodeg(angle(I1))
abs(I1)

%solving for Pin 
Pin = Pout/0.94
PlossT = Pout - Pin %3*Mechanial + 3*core losses.

%finding the output power, calculated.
V2 = Vphase * parr(Xm, X2+R2/s)/RT
I2 = I1 * Xm/(RT)
radtodeg(angle(I2))
abs(I2)

Pag = 3*(I2)^2*R2/s            %air gap power
Pdev = 3*R2/s*(1-s)*(I2)^2    %developed power
Pscl = 3*R1*(I1)^2               %stator copper loss
Prcl = 3*(I2)^2*R2/s            %rotor copper loss
 
Pincalc = 2400/sqrt(3)*I1*cos(angle(I1'))

Pmech = Pincalc - Pout - 3*(I1^2)*R1 -3*R2/s*I2^2
Pcore = Pin - Pout - Pmech -3*(I1^2)*R1 -3*R2/s*I2^2
abs(Pcore)
% finding the Rc
syms rc
s =3.335;
Z1 = R1 + X1;
Z2 = parr(Xm,rc)
Z3 = R2/s + X2;
I1 = Pin/(Vphase); %input into the circuit.

RT = Z1+parr(Z2, Z3)
Vrc = Vphase*parr(parr(rc, Xm), Z3)/RT;

Rc= solve(Vrc^2/RT==Pcore, rc)
