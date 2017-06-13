clc; clear;

%% System parameter

Ts = 1e-4; % Step
Fs = 1/Ts; % Sampling frequency

idref = 10;
iqref = 0;

Udc = 400;
Kpi = 0.1; Kii = 1; % dq control parameters
f1 = 50;

% PLL parameter
xi = 0.707; % Damping coefficient
fPLL = 100; omegaPLL = fPLL * 2 * pi;
V1 = 100;
Ki_PLL = (-2*xi^2*omegaPLL^2+omegaPLL^2*sqrt(1+4*xi^4))/V1;
Kp_PLL = sqrt(4*xi^2*Ki_PLL/V1);

%% Impedance calculation
U1 = 100; % Voltage amp of base frequency
I1 = 10; % Current amp of base frequency
phi1 = 0; % Current angle of base frequency

f = 1:49; % 频率范围

Zp_ideal = (Kpi + Kii./(1j*2*pi*(f-f1)))*Udc/2;
% 理想 PLL 逆变器电阻和电抗
Rp_ideal = real(Zp_ideal); Xp_ideal = imag(Zp_ideal);

HpPLL = (Kp_PLL + Ki_PLL./(1j*2*pi*(f-f1)))*1./(1j*2*pi*(f-f1));
TpPLL = U1*HpPLL./(1+U1*HpPLL);
WpPLL = 1+(Kpi+Kii./(1j*2*pi*(f-f1)))*(Udc/U1)*(I1/2)*exp(1j*phi1);
Zp_Den = 1-1/2 * TpPLL .* WpPLL;
Zp_PLL = Zp_ideal ./ Zp_Den;
% 实际 PLL 逆变器电阻和电抗
Rp_PLL = real(Zp_PLL); Xp_PLL = imag(Zp_PLL);

%% 扫频测量部分

UpScan = 1;
fStart1 = 1;
fEnd1 = 49;
fInt = 1;
fArray1 = fStart1:fInt:fEnd1;
fStart2 = 55;
fEnd2 = 60;
fArray2 = fStart2:fInt:fEnd2;
fArray = fArray1;
% fArray = [fArray1,fArray2];

fNum = length(fArray);

ZpMag = zeros(1,fNum);
ZpPhase = zeros(1,fNum);
ZpR = zeros(1,fNum);
ZpX = zeros(1,fNum);

positiveScan = 1;
if positiveScan == 1
    PhaseB = -120;
    PhaseC = -240;
else
    PhaseB = 120;
    PhaseC = 240;
end

timeSim = 3;

scan_f = true;
if scan_f == true
    % 扫频
    tic;
    for i = 1:fNum % 每个频率扫描一次
        display(strcat('f=',num2str(fArray(i))));
        fpScan = fArray(i);
        
        start_time = 2;
        end_time  = 3;
        
        open('TPI_PLL');
        sim('TPI_PLL');
        
        U = UI.signals(1).values(:,1);
        I = UI.signals(2).values(:,1);
        
        index = start_time*Fs+1:end_time*Fs;
        U = U(index);
        I = I(index);
        dataLen = length(U);
        
        f0 = Fs/dataLen;
        fIndex = fpScan/f0+1;
        
        U_FFT = fft(U);
        
        magU = abs(U_FFT/dataLen);
        magU = magU(1:dataLen/2+1);
        magU(2:end-1) = 2*magU(2:end-1);
        
        phaseU = angle(U_FFT);
        
        I_FFT = fft(I);
        
        magI = abs(I_FFT/dataLen);
        magI = magI(1:dataLen/2+1);
        magI(2:end-1) = 2*magI(2:end-1);
        
        phaseI = angle(I_FFT);
        
        UMag = magU(fIndex);
        UPhase = phaseU(fIndex);
        IMag = magI(fIndex);
        IPhase = phaseI(fIndex) - pi; % 电流相位翻转180度
        
        % 计算阻抗
        if positiveScan == 1
            ZpMag(i) = UMag/IMag;
            ZpPhase(i) =  UPhase - IPhase;
            
            ZpR(i) = ZpMag(i) * cos(ZpPhase(i));
            ZpX(i) = ZpMag(i) * sin(ZpPhase(i));
%         else
%             ZnMag(i) = UMag/IMag;
%             ZnPhase(i) =  UPhase - IPhase;
%             
%             ZnR(i) = ZnMag(i) * cos(ZnPhase(i));
%             ZnX(i) = ZnMag(i) * sin(ZnPhase(i));
        end
    end
    toc;
else
    %% 不扫频
    UpScan = 0; fpScan = 0;
    PhaseB = 0; PhaseC = 0;
end
