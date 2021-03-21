clc
%% user inputs
% L-circuit (shunt Y on the far-receiving side of series Z)
j = sqrt(-1);
k = 1.11; % weird fudge factor for quiver() plots to look right

E2 = 1;
I2 = 0.2;
a = 1.046+j*.0915;  % 1.05 5º
z = .17+j*.98;  % 1 80º
b = 1;
y = (a-1)/z;

%% I2 real

% calculations
E2a = E2*a;
I2z = I2*z;
E1 = E2*a + I2*z;
I1 = E2*y+I2*b;

S2 = E2*conj(I2);
S1 = E1*conj(I1);

% values
P1 = real(S1)
Q1 = imag(S1)
E1mag = abs(E1)

% plots
figure(1)
hold on
quiver(0,0,k*real(S1),k*imag(S1))
quiver(0,0,k*real(E2),k*imag(E2))
quiver(0,0,k*real(I2),k*imag(I2))
quiver(0,0,k*real(E2a),k*imag(E2a))
quiver(real(E2a),imag(E2a),k*real(I2z),k*imag(I2z))
quiver(0,0,k*real(E1),k*imag(E1))
quiver(0,0,k*real(I1),k*imag(I1))

axis([0 1.5 -.5 1])
legend('S1','E2','I2','E2*a','I2z','E1','I1')
hold off

%% I2 changes such that S1 decreases with constant Q
figure(2)
hold on

% plot original I2=real
quiver(0,0,k*real(E2a),k*imag(E2a))
quiver(real(E2a),imag(E2a),k*real(I2z),k*imag(I2z))
quiver(0,0,1.11*real(E1),k*imag(E1))
quiver(0,0,1.11*real(S1),k*imag(S1))


% add .1j to I2
I2 = I2 + .1*j;
E2a = E2*a;
I2z = I2*z;
E1 = E2*a + I2*z;
I1 = E2*y+I2*b;

S2 = E2*conj(I2);
S1 = E1*conj(I1);

% values
P1 = real(S1)/.441
Q1 = imag(S1)/.2131
E1mag = abs(E1)/1.1823

% plot
quiver(0,0,k*real(E2a),k*imag(E2a))
quiver(real(E2a),imag(E2a),k*real(I2z),k*imag(I2z))
quiver(0,0,1.11*real(E1),k*imag(E1))
quiver(0,0,1.11*real(S1),k*imag(S1))

% add .1j to I2
I2 = I2 + .1*j;
E2a = E2*a;
I2z = I2*z;
E1 = E2*a + I2*z;
I1 = E2*y+I2*b;
S1 = E1*conj(I1);

% plot
quiver(0,0,k*real(E2a),k*imag(E2a))
quiver(real(E2a),imag(E2a),k*real(I2z),k*imag(I2z))
quiver(0,0,1.11*real(E1),k*imag(E1))
quiver(0,0,1.11*real(S1),k*imag(S1))

% add .1j to I2
I2 = I2 + .1*j;
E2a = E2*a;
I2z = I2*z;
E1 = E2*a + I2*z;
I1 = E2*y+I2*b;
S1 = E1*conj(I1);

% plot
quiver(0,0,k*real(E2a),k*imag(E2a))
quiver(real(E2a),imag(E2a),k*real(I2z),k*imag(I2z))
quiver(0,0,1.11*real(E1),k*imag(E1))
quiver(0,0,1.11*real(S1),k*imag(S1))


axis([0 1.5 -.5 1])
hold off

