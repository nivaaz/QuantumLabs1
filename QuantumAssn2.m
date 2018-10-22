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
Sy{1} = kron(pauli_x, kron(eye(2), eye(2)));
Sz{1} = kron(pauli_x, kron(eye(2), eye(2)));

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
        h = h + Sz{c}*Sz{d}
    end
end
h = J*h;
[hv, hd] = eig(h)
%%
L = length(h);
l = 1:L;
subplot(2, 1, 1);
plot1 = plot(l, sum(hv), 'o');
xlabel("Energy Number");
ylabel("Eigen Energy");
title("Energy Eigenstate");

subplot(2, 1, 2);
plot2 = plot(l, hd, 'p');
xlabel("Energy Number");
ylabel("Eigen Energy");
title("Eigen Energies");

set([plot2 plot1],'LineWidth',2)

%% Zeeman Interaction
B = 0:1e-6:100e-6;
BS =@(B) 2*pi*gamma_e*B*Sz;

clc
h = 0;
for c = 1:N
    for d = 1:N
        h = h + Sz{c}*Sz{d}
    end
    h = h + BS(B(1))*Sz(c)
end

h = J*h;
[hv, hd] = eig(h)