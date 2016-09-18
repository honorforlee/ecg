load('ecg_noiseBL.mat');

N = length(ecg);
t = [0:N-1]/Fs;

figure(1)
subplot 211
plot(t,ecg);
set(gca,'FontSize',14)
title('ECG with baseline interference','FontSize',16);
xlabel('Time (s)','FontSize',16);
ylabel('ECG signal amplitude','FontSize',16);
xlim([0 4]);

s = 15; 
N1 = Fs * s;
F=([0:N1-1]/N1-0.5)*Fs;
fft_BL_ecg = abs(fft(ecg,N1)).^2;
fft_BL_ecg = fftshift(fft_BL_ecg);

subplot 212
plot(F,fft_BL_ecg);
set(gca,'FontSize',14)
title('ECG with baseline interference spectrum','FontSize',16);
xlabel('Frequency (Hz)','FontSize',16);
ylabel('ECG spectrum amplitude','FontSize',16);
xlim([-30 30]);

% Chebychev filter

[b_cheby1,a_cheby1] = cheby1(1,1,0.01,'high');

[H_cheby,F] = freqz(b_cheby1,a_cheby1,1024,Fs);

figure(2)
subplot 211
plot(F,db(H_cheby)); % Amplitude
set(gca,'FontSize',12)
xlabel('Frequency (Hz)','FontSize',14);
ylabel('H_high filter amplitude |H_{cheby}(f)|','FontSize',14);
title('High-pass filter frequency response','FontSize',14);

subplot 212
plot(F,unwrap(angle(H_cheby)));
set(gca,'FontSize',12)
xlabel('Frequency (Hz)','FontSize',14);
ylabel('H_high filter phase arg(H_{cheby}(f))','FontSize',14); % Phase

ecg_filtered_BL = filter(b_cheby1,a_cheby1,ecg);

figure(3)
subplot 211
plot(t,ecg_filtered_BL);
set(gca,'FontSize',14)
xlabel('Time (s)','FontSize',16); 
ylabel('ECG signal amplitude','FontSize',16);
title('Filtered ECG using cheby1 filter function','FontSize',16);
xlim([0 4]);

s = 15; 
N1 = Fs * s;
F=([0:N1-1]/N1-0.5)*Fs;
fft_BL_filtered_ecg = abs(fft(ecg_filtered_BL,N1)).^2;
fft_BL_filtered_ecg = fftshift(fft_BL_filtered_ecg);

subplot 212
plot(F,fft_BL_filtered_ecg);
set(gca,'FontSize',14)
xlabel('Frequency (Hz)','FontSize',16);
ylabel('ECG spectrum amplitude','FontSize',16);
title('High-pass filtered ECG spectrum','FontSize',16);
xlim([-30 30]);