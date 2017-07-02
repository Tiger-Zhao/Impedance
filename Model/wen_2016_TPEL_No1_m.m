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

Udc = 800;
Us = 690/sqrt(3)*sqrt(2);
L = 0.1 * Lb;
Kp = 0.6; Ki = 15; % PI control paremeters

idref = -sqrt(2) * Ib;
iqref = 0;
I0 = sqrt(idref^2 + iqref^2);

Kp_PLL = 4.46;
Ki_PLL = 991;

Kp_PLL_m = 0.001;
Ki_PLL_m = 0.001;

Z_len = 11;

t_start = 1;
t_end = 2;
f_interval = 1;

f0 = 50;

% Initial value measurement
U_per = Us*0.0;
f_per = 50;
deg_per_A = 0;
deg_per_B = -120;
deg_per_C = 120;
open('wen_2016_TPEL_No1');
sim('wen_2016_TPEL_No1');

Ud0 = Udq0.signals(1).values(end);
Uq0 = Udq0.signals(2).values(end);
Id0 = -Idq0.signals(1).values(end);
Iq0 = -Idq0.signals(2).values(end);
Dd0 = Ddq0.signals(1).values(end);
Dq0 = Ddq0.signals(2).values(end);

U_per = Us*0.01;
%% Impedance calculation
for i = 1:Z_len
    f_dq_i = 200 + f_interval * (i-1);
    
    f_per = f_dq_i + f0;
    deg_per_A = 0;
    deg_per_B = -120;
    deg_per_C = 120;
    
    sim('wen_2016_TPEL_No1');
    
    index = round(t_start*Fs + 1:t_end*Fs);
    n = length(index);
    fs = Fs/n;
    f_index = round(f_dq_i/fs+1);
    
    Vd = Vdq.signals(1).values;
    Vq = Vdq.signals(2).values;
    Id = -Idq.signals(1).values;
    Iq = -Idq.signals(2).values;
    
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
    
    if f_dq_i < f0
        f_per = abs(f_dq_i - f0);
        deg_per_A = 0+90;
        deg_per_B = -120+90;
        deg_per_C = +120+90;
    else
        f_per = f_dq_i - f0;
        deg_per_A = 0+90;
        deg_per_B = +120+90;
        deg_per_C = -120+90;
    end
    
    sim('wen_2016_TPEL_No1');
    
    Vd = Vdq.signals(1).values;
    Vq = Vdq.signals(2).values;
    Id = -Idq.signals(1).values;
    Iq = -Idq.signals(2).values;
    
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
    
    Zdd_dq = [Id1, Iq1;Id2, Iq2]\[Vd1; Vd2];
    Zqd_qq = [Id1, Iq1;Id2, Iq2]\[Vq1; Vq2];
    
    omega_i = 2 * pi * f_dq_i;
    s = 1j * omega_i;
    Z_out = [s*L, -omega_0*L;
        omega_0*L, s*L];
    
    H_PLL = Kp_PLL + Ki_PLL/s;
    G_PLL = H_PLL/(s+Ud0*H_PLL);
    G_PLL_v = [1, Uq0*G_PLL;
        0, 1-Ud0*G_PLL];
    G_PLL_i = [0 Iq0*G_PLL
        0 -Id0*G_PLL];
    G_PLL_d = [0 -Dq0*G_PLL
        0 Dd0*G_PLL];
    G_id = -Udc/(s^2*L^2+omega_0^2*L^2);
    G_id = G_id * [s*L, omega_0*L
        -omega_0*L, s*L];
    G_ci = [Kp+Ki/s, 0
        0, Kp+Ki/s];
    G_ci = -G_ci;
    G_dei = [0, -omega_0*L/Udc;
        omega_0*L/Udc, 0];
    
    Z_cal = (inv(Z_out) + G_id * ((-G_ci+G_dei)...
        * G_PLL_i + G_PLL_d))\...
        (eye(2,2) + G_id * (G_ci-G_dei));
    Z{i,1} = Z_cal;
end