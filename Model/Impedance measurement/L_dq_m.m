clc; clear;

%% System parameter

Ts = 1e-4; % Simulation step
Fs = 1/Ts; % Sampling frequency
f0 = 50; % Base frequency
omega_0 = 2*pi*50; % Base omega

% PLL parameter
Kp_PLL = 0.001;
Ki_PLL = 0.001;

Us = 690/sqrt(3)*sqrt(2);
L = 0.1;

% Impedance length
Z_len = 10;
f_series = zeros(1, Z_len);

% Real and measured impedance result matrix
Z_dd_r = zeros(1, Z_len);
Z_dq_r = zeros(1, Z_len);
Z_qd_r = zeros(1, Z_len);
Z_qq_r = zeros(1, Z_len);
Z_dd_m = zeros(1, Z_len);
Z_dq_m = zeros(1, Z_len);
Z_qd_m = zeros(1, Z_len);
Z_qq_m = zeros(1, Z_len);

% FFT start and end time
t_start = 1;
t_end = 2;

% Impedance scan interval
f_interval = 10;

% Voltage perturbation
U_per = Us * 0.01;
%% Impedance calculation
for i = 1:Z_len
    f_dq_i = 10 + f_interval * (i-1);
    
    % First perturbation
    f_per = f_dq_i + f0;
    deg_per_A = 0;
    deg_per_B = -120;
    deg_per_C = 120;
    
    % Wavedata index
    index = round(t_start*Fs + 1:t_end*Fs);
    n = length(index);
    fs = Fs/n;
    f_index = round(f_dq_i/fs+1);
    
    open('L_dq');
    sim('L_dq');
    
    % Vd, Vq, Id, Iq time series
    Vd = Vdq.signals(1).values;
    Vq = Vdq.signals(2).values;
    Id = Idq.signals(1).values;
    Iq = Idq.signals(2).values;
    
    Vd = Vd(index); Vd = Vd - mean(Vd);
    Vq = Vq(index); Vq = Vq - mean(Vq);
    Vd_FFT = fft(Vd)/n;
    Vq_FFT = fft(Vq)/n;
    Vd1 = Vd_FFT(f_index) * 2;
    Vq1 = Vq_FFT(f_index) * 2;
    
    Id = Id(index); Id = Id - mean(Id);
    Iq = Iq(index); Iq = Iq - mean(Iq);
    Id_FFT = fft(Id)/n;
    Iq_FFT = fft(Iq)/n;
    Id1 = Id_FFT(f_index) * 2;
    Iq1 = Iq_FFT(f_index) * 2;
    
    % Second perturbation
    if f_dq_i < 50
        f_per = abs(f_dq_i - 50);
        deg_per_A = 0+90;
        deg_per_B = -120+90;
        deg_per_C = +120+90;
    else
        f_per = f_dq_i - 50;
        deg_per_A = 0+90;
        deg_per_B = +120+90;
        deg_per_C = -120+90;
    end
    
    sim('L_dq');
    
    Vd = Vdq.signals(1).values;
    Vq = Vdq.signals(2).values;
    Id = Idq.signals(1).values;
    Iq = Idq.signals(2).values;
    
    Vd = Vd(index); Vd = Vd - mean(Vd);
    Vq = Vq(index); Vq = Vq - mean(Vq);
    Vd_FFT = fft(Vd)/n;
    Vq_FFT = fft(Vq)/n;
    Vd2 = Vd_FFT(f_index) * 2;
    Vq2 = Vq_FFT(f_index) * 2;
    
    Id = Id(index); Id = Id - mean(Id);
    Iq = Iq(index); Iq = Iq - mean(Iq);
    Id_FFT = fft(Id)/n;
    Iq_FFT = fft(Iq)/n;
    Id2 = Id_FFT(f_index) * 2;
    Iq2 = Iq_FFT(f_index) * 2;
    
    % dq impedance measurement
    Zdd_dq = [Id1 Iq1;Id2 Iq2]\[Vd1; Vd2];
    Zqd_qq = [Id1 Iq1;Id2 Iq2]\[Vq1; Vq2];
    
    % dq impedance theoretical 
    omega_i = 2 * pi * f_dq_i;
    s = 1j * omega_i;
    Z_out = [s*L, -omega_0*L;
        omega_0*L, s*L];
    
    f_series(i) = f_dq_i;
    Z_dd_r(i) = s*L;
    Z_dq_r(i) = -omega_0*L;
    Z_qd_r(i) = omega_0*L;
    Z_qq_r(i) = s*L;
    
    Z_dd_m(i) = Zdd_dq(1);
    Z_dq_m(i) = Zdd_dq(2);
    Z_qd_m(i) = Zqd_qq(1);
    Z_qq_m(i) = Zqd_qq(2);
end

% Bode plot calculation
Zdd_m_mag = 20*log10(abs(Z_dd_m));
Zdd_r_mag = 20*log10(abs(Z_dd_r));
Zdd_m_angle = angle(Z_dd_m)*180/pi;
Zdd_r_angle = angle(Z_dd_r)*180/pi;

Zdq_m_mag = 20*log10(abs(Z_dq_m));
Zdq_r_mag = 20*log10(abs(Z_dq_r));
Zdq_m_angle = angle(Z_dq_m)*180/pi;
Zdq_r_angle = angle(Z_dq_r)*180/pi;

Zqd_m_mag = 20*log10(abs(Z_qd_m));
Zqd_r_mag = 20*log10(abs(Z_qd_r));
Zqd_m_angle = angle(Z_qd_m)*180/pi;
Zqd_r_angle = angle(Z_qd_r)*180/pi;

Zqq_m_mag = 20*log10(abs(Z_qq_m));
Zqq_r_mag = 20*log10(abs(Z_qq_r));
Zqq_m_angle = angle(Z_qq_m)*180/pi;
Zqq_r_angle = angle(Z_qq_r)*180/pi;

% Plot
subplot(4,2,1);
semilogx(f_series, Zdd_m_mag, f_series, Zdd_r_mag, 'LineWidth', 3);
xlabel('$f/{\textrm {Hz}}$', 'Interpreter', 'LaTeX');
ylabel('$\textrm{Magnitude}/\textrm{dB}$', 'Interpreter', 'LaTeX');
set(gca, 'FontSize', 18);
legend({'$Z_{dd} \, \textrm{measure}$', '$Z_{dd} \, \textrm{model}$'}, 'Interpreter', 'LaTeX');

subplot(4,2,2);
semilogx(f_series, Zdq_m_mag, f_series, Zdq_r_mag, 'LineWidth', 3);
xlabel('$f/{\textrm {Hz}}$', 'Interpreter', 'LaTeX');
ylabel('$\textrm{Magnitude}/\textrm{dB}$', 'Interpreter', 'LaTeX');
legend({'$Z_{dq} \, \textrm{measure}$', '$Z_{dq} \, \textrm{model}$'}, 'Interpreter', 'LaTeX');
set(gca, 'FontSize', 18);

subplot(4,2,3);
semilogx(f_series, Zdd_m_angle, f_series, Zdd_r_angle, 'LineWidth', 3);
xlabel('$f/{\textrm {Hz}}$', 'Interpreter', 'LaTeX');
ylabel('$\textrm{Phase}/\textrm{deg}$', 'Interpreter', 'LaTeX');
legend({'$Z_{dd} \, \textrm{measure}$', '$Z_{dd} \, \textrm{model}$'}, 'Interpreter', 'LaTeX');
set(gca, 'FontSize', 18);

subplot(4,2,4);
semilogx(f_series, Zdq_m_angle, f_series, Zdq_r_angle, 'LineWidth', 3);
xlabel('$f/{\textrm {Hz}}$', 'Interpreter', 'LaTeX');
ylabel('$\textrm{Phase}/\textrm{deg}$', 'Interpreter', 'LaTeX');
legend({'$Z_{dq} \, \textrm{measure}$', '$Z_{dq} \, \textrm{model}$'}, 'Interpreter', 'LaTeX');
set(gca, 'FontSize', 18);

subplot(4,2,5);
semilogx(f_series, Zqd_m_mag, f_series, Zqd_r_mag, 'LineWidth', 3);
xlabel('$f/{\textrm {Hz}}$', 'Interpreter', 'LaTeX');
ylabel('$\textrm{Magnitude}/\textrm{dB}$', 'Interpreter', 'LaTeX');
legend({'$Z_{qd} \, \textrm{measure}$', '$Z_{qd} \, \textrm{model}$'}, 'Interpreter', 'LaTeX');
set(gca, 'FontSize', 18);

subplot(4,2,6);
semilogx(f_series, Zqd_m_mag, f_series, Zqd_r_mag, 'LineWidth', 3);
xlabel('$f/{\textrm {Hz}}$', 'Interpreter', 'LaTeX');
ylabel('$\textrm{Magnitude}/\textrm{dB}$', 'Interpreter', 'LaTeX');
legend({'$Z_{qq} \, \textrm{measure}$', '$Z_{qq} \, \textrm{model}$'}, 'Interpreter', 'LaTeX');
set(gca, 'FontSize', 18);

subplot(4,2,7);
semilogx(f_series, Zqq_m_angle, f_series, Zqq_r_angle, 'LineWidth', 3);
xlabel('$f/{\textrm {Hz}}$', 'Interpreter', 'LaTeX');
ylabel('$\textrm{Phase}/\textrm{deg}$', 'Interpreter', 'LaTeX');
legend({'$Z_{qd} \, \textrm{measure}$', '$Z_{qd} \, \textrm{model}$'}, 'Interpreter', 'LaTeX');
set(gca, 'FontSize', 18);

subplot(4,2,8);
semilogx(f_series, Zqq_m_angle, f_series, Zqq_r_angle, 'LineWidth', 3);
xlabel('$f/{\textrm {Hz}}$', 'Interpreter', 'LaTeX');
ylabel('$\textrm{Phase}/\textrm{deg}$', 'Interpreter', 'LaTeX');
legend({'$Z_{qq} \, \textrm{measure}$', '$Z_{qq} \, \textrm{model}$'}, 'Interpreter', 'LaTeX');
set(gca, 'FontSize', 18);