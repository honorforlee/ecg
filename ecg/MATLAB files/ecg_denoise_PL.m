load('ecg_noisePL.mat');

N = length(ecg);
t = [0:N-1]/Fs;

figure(1)
subplot 211
plot(t,ecg);
set(gca,'FontSize',14)
title('ECG with powerline interference','FontSize',16);
xlabel('Time (s)','FontSize',16);
ylabel('ECG signal amplitude','FontSize',16);
xlim([0 4]);

s = 15; 
N1 = Fs * s;
F=([0:N1-1]/N1-0.5)*Fs;
fft_PL_ecg = abs(fft(ecg,N1)).^2;
fft_PL_ecg = fftshift(fft_PL_ecg);

subplot 212
plot(F,fft_PL_ecg);
set(gca,'FontSize',14)
title('ECG with powerline interference spectrum','FontSize',16);
xlabel('Frequency (Hz)','FontSize',16);
ylabel('ECG spectrum amplitude','FontSize',16);
xlim([-30 30]);

% Fir1 filter

[b_fir1,a_fir1] = fir1(50,0.01,'low');

[H_fir,F] = freqz(b_fir1,a_fir1,1024,Fs);

figure(2)
subplot 211
plot(F,db(H_fir)); % Amplitude
set(gca,'FontSize',12)
xlabel('Frequency (Hz)','FontSize',14);
ylabel('H_low filter amplitude |H_{fir1}(f)|','FontSize',14);
title('Low-pass filter frequency response','FontSize',14);

subplot 212
plot(F,unwrap(angle(H_fir)));
set(gca,'FontSize',12)
xlabel('Frequency (Hz)','FontSize',14);
ylabel('H_low filter phase arg(H_{fir1}(f))','FontSize',14); % Phase

ecg_filtered_PL = filter(b_fir1,a_fir1,ecg);

figure(3)
subplot 211
plot(t,ecg_filtered_PL);
set(gca,'FontSize',14)
xlabel('Time (s)','FontSize',16); 
ylabel('ECG signal amplitude','FontSize',16);
title('Filtered ECG using fir1 filter function','FontSize',16);
xlim([0 4]);

s = 15; 
N1 = Fs * s;
F=([0:N1-1]/N1-0.5)*Fs;
fft_PL_filtered_ecg = abs(fft(ecg_filtered_PL,N1)).^2;
fft_PL_filtered_ecg = fftshift(fft_PL_filtered_ecg);

subplot 212
plot(F,fft_PL_filtered_ecg);
set(gca,'FontSize',14)
xlabel('Frequency (Hz)','FontSize',16);
ylabel('ECG spectrum amplitude','FontSize',16);
title('Low-pass filtered ECG spectrum','FontSize',16);
xlim([-30 30]);
