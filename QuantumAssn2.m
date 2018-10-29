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
     if d~=c %when they're not equal
        h = h + Sz{c}*Sz{d} %find new hamiltonian
     end
    end
end
h = J*h;
[hv, hd] = eig(h)

%%
L = length(h);
l = 1:L;
subplot(2, 1, 1);
plot1 = plot(l, hv, 'o');
xlabel("Energy Number");
ylabel("Eigen Energy");
title("Energy Eigenstate");

subplot(2, 1, 2);
plot2 = plot(l, hd,'o');
xlabel("Energy Number");
ylabel("Eigen Energy");
title("Eigen Energies");

set([plot2 plot1],'LineWidth',2)

%% Zeeman Interaction
clc
Sztot = Sz{1} + Sz{2} + Sz{3};

Bz = 0:1e-7:100e-6;
BS =@(B)2*pi*gamma_e*B*Sztot;
len = length(Bz);

Engv = zeros(8, 8*len);
Engd = zeros(8, 8*len);

for c =1:len
    H = 0;
    H = h + BS(Bz(c));% hamiltonian, at Bz
    [engv, engd] = eig(H);  %getting eig of new hamiltonian
    if(c==1)
       [Engv, Engd] = eig(H);    %initialalising
    else
        Engv = horzcat(Engv, engv);     %appending new hamiltonian
        Engd = horzcat(Engd, engd);
    end
end

% plot(Bplot, Engd, 'o');
% title('Zeeman Energy');
% xlabel('Field Strength (T)');
% ylabel('Energy (eV)?');

%% PLOTS WITH ZEEMAN INTERACITON
close all;

Bplot = 0
for c = 1:len-1
    if c==1
        Bplot(1:8) = Bz(c);
    else
        Bplot(c*8:8*(c+1)) = Bz(c);
    end
end

figure;
title("Zeeman Interaction")
subplot(2, 1, 1);
plot1 = plot(Bplot, Engv, 'o');
xlabel("Bz");ylabel("Eigen Energy");title("Energy Eigenstate");

subplot(2, 1, 2);
plot2 = plot(Bplot, Engd, 'o');
xlabel("Bz");ylabel("Eigen Energy");title("Eigen Energies");

set([plot2 plot1],'LineWidth',1);

%% 1.4 solenoid
clc
d = 10;
Bz = 0:1e-7:100e-6;

Bxyz = @(Bo, x, y, z) Bo * (1 - (x+y)/d);
BS =@(B, x, y, z) 2*pi*gamma_e*Bxyz(x, y, z)*Sztot;
len = length(Bz);
Engv0 = 0;Engd0 = 0;Engv1 = 0;Engd1 = 0;Engv2 = 0;Engd2 = 0;
 
% 0 0 
a = 0, b = 0, c = 0;
for c =1:len
    H = 0;
    H = h + BS(Bz(c), a, b, c);
    if(c==1)
       [Engv0, Engd0] = eig(H);
    else
        [engv, engd] = eig(H);
        Engv0 = horzcat(Engv0, engv);
        Engd0 = horzcat(Engd0, engd);
    end
end

%1 0 
a = 1, b = 0, c = 0;
for c =1:len
    H = 0;
    H = h + BS(Bz(c), a, b, c);    
    if(c==1)
       [Engv1, Engd1] = eig(H);
    else
        [engv, engd] = eig(H);
        Engv1 = horzcat(Engv1, engv);
        Engd1 = horzcat(Engd1, engd);
    end
end

%0.5, sqrt(3/4)
a = 0.5, b = sqrt(3)/2, c = 0;
for c =1:len
    H = 0;
    H = h + BS(Bz(c), a, b, c)
    
    if(c==1)
       [Engv2, Engd2] = eig(H);
    else
        [engv, engd] = eig(H);
        Engv2 = horzcat(Engv2, engv);
        Engd2 = horzcat(Engd2, engd);
    end
end

Bplot = 0
for c = 1:len-1
    if c==1
        Bplot(1:8) = Bz(c);
    else
        Bplot(c*8:8*(c+1)) = Bz(c);
    end
end

subplot(3, 2, 1);
plot(Bplot, Engv0, 'o');
ylabel("Eigen-Vector");title("Spin at  0,0");

subplot(3, 2, 2);
plot(Bplot, Engd0, 'o');
ylabel("Eigen-Energy");

subplot(3, 2, 3);
plot(Bplot, Engv1, 'o');
ylabel("Eigen-Vector");title("Spin at  1,0");

subplot(3, 2, 4);
plot(Bplot, Engd1, 'o');
ylabel("Eigen-Energy");

subplot(3, 2, 5);
plot(Bplot, Engv2, 'o');
ylabel("Eigen-Vector");title("Spin at  1/2,sqrt(3/4)");

subplot(3, 2, 6);
plot(Bplot, Engd2, 'o');
ylabel("Eigen-Energy");

title('Field with Distance');
%% Heisenberg form 
clc
for c = 1:N     %setting up the hamiltonian
    for d = 1:N
    if (c~=d)
        h = h + Sx{d}*Sx{c} + Sy{d}*Sy{c} + Sz{d}*Sz{c}
    end
    end
end
h = j*h; %multiplying by antiferro

% Run for Zeeman splittng , WITH FIELD 
Engv = zeros(8, 8*len);
Engd = zeros(8, 8*len);

Bz = 0:1e-7:100e-6;
BS =@(B)2*pi*gamma_e*B*Sztot;
len = length(Bz);

for c =1:len        %   FOR NO FIELD FORCED.
    H = 0;
    H = h + BS(Bz(c));              % hamiltonian, at Bz
    [engv, engd] = eig(H);          %getting eig of new hamiltonian
    if(c==1)
       [Engv, Engd] = eig(H);       %initialalising
    else
        Engv = horzcat(Engv, engv); %appending new hamiltonian
        Engd = horzcat(Engd, engd);
    end
end

figure;
title("Zeeman Interaction");
subplot(2, 2, 1);
plot1 = plot(Bplot, Engv, 'o');
xlabel("Bz");ylabel("Zeeman Energy");title("Energy Eigenstate");

subplot(2, 2, 2);
plot2 = plot(Bplot, Engd, 'o');
xlabel("Bz");ylabel("Zeeman Energy");title("Eigen Energies");
set([plot2 plot1],'LineWidth',1);

[Ev, Ed] = eig(h);
subplot(2, 2, 3);
plot1 = plot(1:length(Ev), Ev, 'o');
xlabel("Bz");ylabel("Ising Energy");title("Energy Eigenstate");

subplot(2, 2, 4);
plot2 = plot(1:length(Ed), Ed, 'o');
xlabel("Bz");ylabel("Ising Energy");title("Eigen Energies");

%%1.6
