clear 
clc
hbar = 1.0545718e-34;
J = 1e6;
N = 3; %number of spins
gamma_e = 29e9;

pauli_x = hbar/2*[0 1; 1 0]
pauli_y = hbar/2*[0 -j; j 0]
pauli_z = hbar/2*[1 0;0 -1]

Sx = {}; Sy = {}; Sz = {};

% for c = 1:N
%     a = eye(2);
%     for d = 1:N-1
%         if (c==d)
%             a = kron(pauli_x, a)
%         else
%            a = kron(eye(2), a) 
%         end
%     end
%     S{c} = a
% end
% %%
% clc
% kron(eye(2), kron(pauli_x, eye(2)))
% S{2}
%% SPIN OPERATORS.
clc
Sx{1} = kron(pauli_x, kron(eye(2), eye(2)));
Sy{1} = kron(pauli_y, kron(eye(2), eye(2)));
Sz{1} = kron(pauli_z, kron(eye(2), eye(2)));

Sx{2} = kron(eye(2), kron(pauli_x, eye(2)));
Sy{2} = kron(eye(2), kron(pauli_y, eye(2)));
Sz{2} = kron(eye(2), kron(pauli_z, eye(2)));

Sx{3} = kron(eye(2), kron(eye(2), pauli_x));
Sy{3} = kron(eye(2), kron(eye(2), pauli_y));
Sz{3} = kron(eye(2), kron(eye(2), pauli_z));

%% Ising form 
clc
h = 0;
for c = 1:N
    for d = 1:N
        if d==c
        h = h + Sz{c}*Sz{d}
        end
    end
end
h = J*h;
[hv, hd] = eig(h)

%%
L = length(h);
l = 1:L;
subplot(2, 1, 1);
plot1 = plot(l, hv, 'p');
xlabel("Energy Number");
ylabel("Eigen Energy");
title("Energy Eigenstate");

subplot(2, 1, 2);
plot2 = plot(l, hd, 'o');
xlabel("Energy Number");
ylabel("Eigen Energy");
title("Eigen Energies");

set([plot2 plot1],'LineWidth',2)

%% Zeeman Interaction
clc

Sztot = Sz{1} + Sz{2} + Sz{3};

Bz = 0:1e-9:100e-6;

BS =@(B) 2*pi*gamma_e*B*Sztot;

len = length(B);

for c =1:len
    H = 0;
    H = h + BS(B(c))
    [engv, engd] = eig(H);
    if(c==1)
       [Engv, Engd] = eig(H)
    else
        Engv = horzcat(Engv, engv);
        Engd = horzcat(Engd, engd);
    end
end

%% PLOTS WITH ZEEMAN INTERACITON
close all;
L = length(Engv);
l = 1:L;

figure;
title("Zeeman Interation")
subplot(2, 1, 1);
plot1 = plot(l, Engv, 'o');
xlabel("Energy Number");ylabel("Eigen Energy");title("Energy Eigenstate");

subplot(2, 1, 2);
plot2 = plot(l, Engd, 'o');
xlabel("Energy Number");ylabel("Eigen Energy");title("Eigen Energies");

set([plot2 plot1],'LineWidth',1);

%% 1.4 solenoid
clc
d = 10;
Bz = 0:1e-9:100e-6;

t = linspace(0, 100e-6,808); %for plotting, end val will change if step changes!

Bxyz = @(Bo, x, y, z) Bo * (1 - (x+y)/d);

BS =@(B, x, y, z) 2*pi*gamma_e*Bxyz(x, y, z)*Sztot;

len = length(B);

% 0 0 
a = 0, b = 0, c = 0;
for c =1:len
    H = 0;
    H = h + BS(B(c), a, b, c)
    
    if(c==1)
       [Engv0, Engd0] = eig(H);
    else
        [engv, engd] = eig(H);
        Engv0 = horzcat(Engv0, engv);
        Engd0 = horzcat(Engd0, engd);
    end
end

subplot(3, 2, 1);
plot(t, Engv0, 'o');
ylabel("Eigen-Vector");title("Spin at  0,0");

subplot(3, 2, 2);
plot(t, Engd0, 'o');
ylabel("Eigen-Energy");


%1 0 
a = 1, b = 0, c = 0;
for c =1:len
    H = 0;
    H = h + BS(B(c), a, b, c)
    
    if(c==1)
       [Engv1, Engd1] = eig(H);
    else
        [engv, engd] = eig(H);
        Engv1 = horzcat(Engv1, engv);
        Engd1 = horzcat(Engd1, engd);
    end
end

subplot(3, 2, 3);
plot(t, Engv1, 'o');
ylabel("Eigen-Vector");title("Spin at  1,0");

subplot(3, 2, 4);
plot(t, Engd1, 'o');
ylabel("Eigen-Energy");

%0.5, sqrt(3/4)
a = 0.5, b = sqrt(3)/2, c = 0;
for c =1:len
    H = 0;
    H = h + BS(B(c), a, b, c)
    
    if(c==1)
       [Engv2, Engd2] = eig(H);
    else
        [engv, engd] = eig(H);
        Engv2 = horzcat(Engv2, engv);
        Engd2 = horzcat(Engd2, engd);
    end
end
subplot(3, 2, 5);
plot(t, Engv2, 'o');
ylabel("Eigen-Vector");title("Spin at  1/2,sqrt(3/4)");
subplot(3, 2, 6);
plot(t, Engd2, 'o');
ylabel("Eigen-Energy");

%% Heisenberg form 
for c = 1:N
    for d = 1:N
    if (c==d)
        h = h + Sx{d}*Sx{c} + Sy{d}*Sy{c} + Sz{d}*Sz{c}
    end
end
