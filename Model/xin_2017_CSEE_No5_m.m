clc; clear;

%% System parameter

Ts = 1e-4; % Simulation step
Fs = 1/Ts; % Sampling frequency
omega_0 = 2*pi*50; % Base omega

Sb = 500e3; % Base VA
Ub = 0.69e3; % Base V
Ib = Sb/(sqrt(3)*Ub); % Base I
Zb = Ub^2/Sb; % Base Z
Lb = Zb/omega_0; % Base L

Udc = 1;
Us = 690/sqrt(3)*sqrt(2);
U0 = Us * 0.9837;
Lf = 0.2 * Lb;
L_line = 0.1 * Lb;
Kp2 = 0.6; Ki2 = 15; % PI control paremeters
f1 = 50;

idref = sqrt(2) * Ib;
iqref = 0;
I0 = sqrt(idref^2 + iqref^2);

Kp_PLL = 2.5;
Ki_PLL = 3020;
T_FF = 0.001;

% Nyquist of -sYg4(s)
tf1 = -tf([Kp2, Ki2], [Lf, Kp2, Ki2]);
tf2 = tf(I0*[Kp_PLL, Ki_PLL], [1, U0*Kp_PLL, U0*Ki_PLL]);
tf3 = -tf([1, 0], 1);
tf = tf1 * tf2 * tf3; 
tf = tf * Lb; % pu(L) Nyquist

% nyquist(tf);
% 
% omega = 1:1:1500;
% omega = omega * 1j * 2 * pi;
% num_val = polyval(tf.num{1},omega);
% den_val = polyval(tf.den{1},omega);
% result = num_val./den_val;
% [Y,I] = min(abs(imag(result)));