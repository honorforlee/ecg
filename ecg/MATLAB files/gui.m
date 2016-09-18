function varargout = gui(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

function gui_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);

function varargout = gui_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function edit_time_min_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Load ECG signal
function load_ecg_file_Callback(hObject, eventdata, handles)
 
[filename,pathname] = uigetfile({'*.mat'},'Load ECG signal');
path = [pathname filename];

if (filename == 0)
   disp('Cancel')
else
    disp(['ECG loaded: ', fullfile(pathname,filename)])
    Struct = load(filename);
    handles.filename = filename;
end
    
 guidata(hObject,handles)

 
% Display ECG signal in time interval
function ecg_time_display_Callback(hObject, eventdata, handles)

filename=handles.filename;
Struct = load(filename);
Fs = Struct.Fs;
ecg = Struct.ecg;
t = [0:length(ecg)-1]/Fs;

time_min = str2double(get(handles.edit_time_min,'String'));
time_max = str2double(get(handles.edit_time_max,'String'));
time_interval = time_min:1/Fs:time_max;
ech_interval = (1+time_min)*Fs:(1+time_max)*Fs;

if isnan(time_min) | isnan(time_max) % Default behavior if no input
   plot(handles.axes1,t,ecg);
   xlim([0,5])
else
    plot(handles.axes1,time_interval,ecg(ech_interval));
end

title(handles.axes1,'ECG Signal','FontName','Helvetica','FontSize',16);
xlabel(handles.axes1,'Time (s)','FontName','Helvetica','FontSize',16);
ylabel(handles.axes1,'ECG Signal Amplitude','FontName','Helvetica','FontSize',16);

if isnan(time_min) | isnan(time_max)
    ech_interval = Fs:6*Fs;
    threshold = (max(ecg([1:5000]))-mean(ecg([1:5000])))/3+mean(ecg);
    [peak,peak_loc] = findpeaks(ecg(ech_interval),'MinPeakHeight',threshold);
else
    threshold = (max(ecg([1:5000]))-mean(ecg([1:5000])))/3+mean(ecg);
    [peak,peak_loc] = findpeaks(ecg(ech_interval),'MinPeakHeight',threshold);
end

if ~isempty(peak)
    
    for i=1:length(peak_loc)-1
        d(i) = (peak_loc(i+1)-peak_loc(i))/Fs;
        d_mean = mean(d(i));
    end

    for i=1:length(d)

        bpm(i)=60/d(i);
        bpm=round(mean(bpm(i)));

        if (bpm > 800)
            state = 'Warning ! Heart rate too high';
            color = 'r';
        elseif (bpm < 20)
            state = 'Warning ! Heart rate too low';
            color = 'r';
        elseif (bpm < 600) & (bpm > 240)
            state = 'Ventricular Fibrillation';
            color = 'r';
        elseif (bpm > 110)
            state = 'Tachycardia';
            color = 'r';
        elseif (bpm < 60);
            state = 'Bradycardia';
            color = 'r';
        else
            state = 'Normal';
            color = [0,0.498,0];
        end
    end
    
    set(handles.bpm_display,'String',bpm,'ForegroundColor',color);
    set(handles.anomaly_display,'String',state,'ForegroundColor',color);
    
end

handles.d = d;
handles.peak_loc = peak_loc;
handles.t = t;
handles.time_min = time_min;
handles.time_max = time_max;
handles.ech_interval = ech_interval;
handles.time_interval = time_interval;

guidata(hObject,handles) 

% Choose Min Time input
function edit_time_min_Callback(hObject, eventdata, handles)

% Choose Max Time input
function edit_time_max_Callback(hObject, eventdata, handles)

% Display ECG spectrum in frequency interval
function ecg_spectrum_display_Callback(hObject, eventdata, handles)

filename=handles.filename;
Struct = load(filename);
Fs = Struct.Fs;
ecg = Struct.ecg;

s = 5;
N = Fs * s;
F=([0:N-1]/N-0.5)*Fs;
fft_ecg = abs(fft(ecg([1:N]),N)).^2;
fft_ecg = fftshift(fft_ecg);

freq_min = str2double(get(handles.edit_freq_min,'String'));
freq_max = str2double(get(handles.edit_freq_max,'String'));
freq_interval = freq_min:freq_max;

if isnan(handles.time_min) | isnan(handles.time_max) % Default behavior if no input
   plot(handles.axes2,F,fft_ecg);
   xlim(handles.axes2,[-50 50]);
    
   fft_threshold = mean(fft_ecg);
   peak_search_start = floor(length(fft_ecg)/2)+5;
   peak_search_end = floor(length(fft_ecg)/2)+5+100;
   [fft_peak fft_peak_loc] = findpeaks(fft_ecg([peak_search_start:peak_search_end]),'MinPeakHeight',fft_threshold);
   
   if isempty(fft_peak)
      set(handles.bpm_fft_display,'String','Error');
   end
   
   fft_peak_max = fft_peak(1);
   first_fft_peak_ech_loc = ceil(length(fft_ecg)/2 + fft_peak_loc(1));
   first_fft_peak_loc = F(first_fft_peak_ech_loc+3);
    
elseif isnan(freq_min) | isnan(freq_max)
       N1 = length(handles.ech_interval);
       F1=([0:N1-1]/N1-0.5)*Fs;
       fft_ecg = abs(fft(ecg(handles.ech_interval),N1)).^2;
       fft_ecg = fftshift(fft_ecg);
       plot(handles.axes2,F1,fft_ecg);
       xlim(handles.axes2,[-50 50])

       fft_threshold = mean(fft_ecg);
       peak_search_start = floor(length(fft_ecg)/2)+5;
       peak_search_end = floor(length(fft_ecg)/2)+5+100;
       [fft_peak fft_peak_loc] = findpeaks(fft_ecg([peak_search_start:peak_search_end]),'MinPeakHeight',fft_threshold);
       
       if isempty(fft_peak)
          set(handles.bpm_fft_display,'String','Error');
       end
    
       fft_peak_max = fft_peak(1);
       first_fft_peak_ech_loc = ceil(length(fft_ecg)/2 + fft_peak_loc(1));
       first_fft_peak_loc = F1(first_fft_peak_ech_loc+3);
  
else
       N1 = length(handles.ech_interval);
       F1=([0:N1-1]/N1-0.5)*Fs;
       fft_ecg = abs(fft(ecg(handles.ech_interval),N1)).^2;
       fft_ecg = fftshift(fft_ecg);
       plot(handles.axes2,F1,fft_ecg);
       xlim(handles.axes2,[freq_min freq_max])

       fft_threshold = mean(fft_ecg);
       peak_search_start = floor(length(fft_ecg)/2)+5;
       peak_search_end = floor(length(fft_ecg)/2)+5+100;
       [fft_peak fft_peak_loc] = findpeaks(fft_ecg([peak_search_start:peak_search_end]),'MinPeakHeight',fft_threshold);
       
       if isempty(fft_peak)
          set(handles.bpm_fft_display,'String','Error');
       end
    
       fft_peak_max = fft_peak(1);
       first_fft_peak_ech_loc = ceil(length(fft_ecg)/2 + fft_peak_loc(1));
       first_fft_peak_loc = F1(first_fft_peak_ech_loc+3);
    
end

bpm_fft = round(first_fft_peak_loc*60);

if (bpm_fft > 800)
    state = 'Warning ! Heart rate too high';
    color = 'r';
elseif (bpm_fft < 20)
    state = 'Warning ! Heart rate too low';
    color = 'r';
elseif (bpm_fft < 600) & (bpm_fft > 240)
    state = 'Ventricular Fibrillation';
    color = 'r';
elseif (bpm_fft > 110)
    state = 'Tachycardia';
    color = 'r';
elseif (bpm_fft < 60)
    state = 'Bradycardia';
    color = 'r';
else
    state = 'Normal';
    color = [0,0.498,0];
end

set(handles.bpm_fft_display,'String',bpm_fft,'ForegroundColor',color);

title(handles.axes2,'ECG Spectrum','FontName','Helvetica','FontSize',16);
xlabel(handles.axes2,'Frequency (Hz)','FontName','Helvetica','FontSize',16);
ylabel(handles.axes2,'ECG Spectrum Amplitude','FontName','Helvetica','FontSize',16);

guidata(hObject,handles) 

% Choose Min Freq input
function edit_freq_min_Callback(hObject, eventdata, handles)

% Choose Max Freq input
function edit_freq_max_Callback(hObject, eventdata, handles)

% PQRST waves detection
function waves_detect_Callback(hObject, eventdata, handles)

filename=handles.filename;
Struct = load(filename);
Fs = Struct.Fs;
ecg = Struct.ecg;

N = length(ecg);
t = [0:N-1]/Fs;

%% Low-pass filter
l_b = [1 zeros(1,5) -2 zeros(1,5) 1];
l_a = [1,-2,1];

ecg_low_filtered = filter(l_b,l_a,ecg); % ECG signal after low-pass filtering

%% High-pass filter
h_b = [-1 zeros(1,15) 32,-32 zeros(1,14) 1];
h_a = [1,1];

ecg_high_filtered = filter(h_b,h_a,ecg_low_filtered); % ECG signal after high-pass filtering
ecg_filtered = ecg_high_filtered; % ECG signal after band-pass filtering

%% Deriving

% Differentiation filter
b_diff = [Fs/8,2*Fs/8,0,-2*Fs/8,-Fs/8];
a_diff = [1];

ecg_derived = filter(b_diff,a_diff,ecg_filtered); % ECG signal derived

%% Squaring

ecg_squared = abs(ecg_derived).^2;

%% Moving-window integration

Win_width = 40;

% Moving window integration filter
b_MWI = ones(1,Win_width);
a_MWI = [Win_width];

ecg_MWI = filter(b_MWI,a_MWI,ecg_squared);

%% Thresholding

ecg_MWI_deriv = diff(ecg_MWI);
t2 = [0:N-2]/Fs;

max_vect = [max(abs(ecg_MWI_deriv([1*Fs:10*Fs]))) max(abs(ecg_MWI_deriv([40*Fs:50*Fs]))) max(abs(ecg_MWI_deriv([90*Fs:100*Fs]))) ];
threshold = mean(max_vect)/1.8;

[peak_max,peak_loc_max] = findpeaks(ecg_MWI_deriv,'MinPeakHeight',threshold);
ecg_MWI_deriv_inv = -ecg_MWI_deriv;
[peak_min,peak_loc_min] = findpeaks(ecg_MWI_deriv_inv,'MinPeakHeight',threshold);

ecg_MWI_slope_begin = peak_loc_max;
ecg_MWI_slope_end = peak_loc_min;

for i = 1:min(length(ecg_MWI_slope_begin),length(ecg_MWI_slope_end))
    for j=ecg_MWI_slope_begin(i):ecg_MWI_slope_end(i)
        if (ecg_MWI(j+1) > ecg_MWI(j))
            R_peak_loc_range_mean(j) = (ecg_MWI_slope_begin(i)+ecg_MWI_slope_end(i))/2;
        end
    end
end

R_peak_loc_range = R_peak_loc_range_mean(find(R_peak_loc_range_mean > 0));

for i = 1:length(R_peak_loc_range)-1
    if R_peak_loc_range(i+1)-R_peak_loc_range(i)>50
        R_peak_loc(i) = R_peak_loc_range(i);
    end
end

R_peak_loc = round(R_peak_loc(find(R_peak_loc>0)));

R_peak_loc_delay = R_peak_loc - Win_width - 2;

for i=1:length(R_peak_loc_delay)
    if (ecg(R_peak_loc_delay(i)) < 0)
        R_peak_loc_delay(i) = R_peak_loc_delay(i) - 8;
    end
end

%% Q and S waves detection

for i=1:length(R_peak_loc_delay)-1
    [S_peak(i) S_peak_loc(i)] = min(ecg([R_peak_loc_delay(i):R_peak_loc_delay(i)+15]));
    S_peak_loc(i)=S_peak_loc(i)+R_peak_loc_delay(i);
end

for i=1:length(R_peak_loc_delay)-1
    [Q_peak(i) Q_peak_loc(i)] = min(ecg([R_peak_loc_delay(i)-15:R_peak_loc_delay(i)]));
    Q_peak_loc(i)=R_peak_loc_delay(i)-15+Q_peak_loc(i);
end

%% P and T waves detection

for i=1:length(R_peak_loc_delay)-1
    [T_peak(i) T_peak_loc(i)] = max(ecg([R_peak_loc_delay(i)+30:R_peak_loc_delay(i)+0.7*(R_peak_loc_delay(i+1)-R_peak_loc_delay(i))]));
    T_peak_loc(i)=T_peak_loc(i)+R_peak_loc_delay(i)+30;
end

for i=1:length(R_peak_loc_delay)-1
    Range(i) = round(R_peak_loc_delay(i)+0.7*(R_peak_loc_delay(i+1)-R_peak_loc_delay(i)));
end

for i=1:length(R_peak_loc_delay)-1
    [P_peak(i) P_peak_loc(i)] = max(ecg([Range(i):R_peak_loc_delay(i+1)-20]));
    P_peak_loc(i) = P_peak_loc(i)+Range(i);
end

if isnan(handles.time_min) | isnan(handles.time_max) % Default behavior if no input
   plot(handles.axes1,t,ecg);
   xlim([0,5])
   
   hold on;
   plot(P_peak_loc/Fs,ecg(P_peak_loc),'cv','MarkerFaceColor','c','MarkerSize',8);
   plot(Q_peak_loc/Fs,ecg(Q_peak_loc),'gv','MarkerFaceColor','g','MarkerSize',8);
   plot(R_peak_loc_delay/Fs,ecg(R_peak_loc_delay),'rv','MarkerFaceColor','r','MarkerSize',8);
   plot(S_peak_loc/Fs,ecg(S_peak_loc),'bv','MarkerFaceColor','b','MarkerSize',8);
   plot(T_peak_loc/Fs,ecg(T_peak_loc),'mv','MarkerFaceColor','m','MarkerSize',8);
   legend('ECG signal','P waves','Q waves','R waves','S waves','T waves');
   hold off;
else
    plot(handles.axes1,handles.time_interval,ecg(handles.ech_interval));
    xlim([handles.time_min, handles.time_max])
   
    hold on;
    plot(P_peak_loc/Fs-1,ecg(P_peak_loc),'cv','MarkerFaceColor','c','MarkerSize',8);
    plot(Q_peak_loc/Fs-1,ecg(Q_peak_loc),'gv','MarkerFaceColor','g','MarkerSize',8);
    plot(R_peak_loc_delay/Fs-1,ecg(R_peak_loc_delay),'rv','MarkerFaceColor','r','MarkerSize',8);
    plot(S_peak_loc/Fs-1,ecg(S_peak_loc),'bv','MarkerFaceColor','b','MarkerSize',8);
    plot(T_peak_loc/Fs-1,ecg(T_peak_loc),'mv','MarkerFaceColor','m','MarkerSize',8);
    legend('ECG signal','P waves','Q waves','R waves','S waves','T waves');
    hold off;
end

handles.R_peak_loc_delay = R_peak_loc_delay;
   
title(handles.axes1,'ECG Signal','FontName','Helvetica','FontSize',16);
xlabel(handles.axes1,'Time (s)','FontName','Helvetica','FontSize',16);
ylabel(handles.axes1,'ECG Signal Amplitude','FontName','Helvetica','FontSize',16);

guidata(hObject,handles) 

% Baseline & Powerline denoising
function ecg_denoise_BL_PL_Callback(hObject, eventdata, handles)

filename=handles.filename;
Struct = load(filename);
Fs = Struct.Fs;
ecg = Struct.ecg;

N = length(ecg);
t = [0:N-1]/Fs;

[b_cheby1,a_cheby1] = cheby1(1,1,0.01,'high');
ecg_filtered_BL = filter(b_cheby1,a_cheby1,ecg);

[b_fir1,a_fir1] = fir1(50,0.01,'low');
ecg_filtered_BL_PL = filter(b_fir1,a_fir1,ecg_filtered_BL);

if isnan(handles.time_min) | isnan(handles.time_max) % Default behavior if no input
    plot(handles.axes1,t,ecg_filtered_BL_PL);
    xlim([0,5])
else
    plot(handles.axes1,handles.time_interval,ecg_filtered_BL_PL(handles.ech_interval));
end

title(handles.axes1,'ECG Signal','FontName','Helvetica','FontSize',16);
xlabel(handles.axes1,'Time (s)','FontName','Helvetica','FontSize',16);
ylabel(handles.axes1,'ECG Signal Amplitude','FontName','Helvetica','FontSize',16);

guidata(hObject,handles) 

% Ectopic beats detection
function ectopic_detect_Callback(hObject, eventdata, handles)
%% Ectopic beats
filename=handles.filename;
Struct = load(filename);
Fs = Struct.Fs;
ecg = Struct.ecg;

for i=1:length(handles.R_peak_loc_delay)-1
    distance_R(i) = (handles.R_peak_loc_delay(i+1)-handles.R_peak_loc_delay(i))/Fs;
end

time_threshold = 0.5;

for i=1:length(handles.R_peak_loc_delay)-2
    if (distance_R(i+1) - distance_R(i)) > time_threshold
       ectopic_beats(i) = handles.R_peak_loc_delay(i+1)/Fs;
    end
end
    
ectopic_beats = ectopic_beats(find(ectopic_beats > 0));

if isempty(ectopic_beats)
   set(handles.ectopic_display,'String','No Ectopic beats detected');
else
   set(handles.ectopic_display,'String','Ectopic beats detected');
end

guidata(hObject,handles) 

function edit_freq_min_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function figure1_ResizeFcn(hObject, eventdata, handles)

function edit_time_max_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_freq_max_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function bpm_display_Callback(hObject, eventdata, handles)

function bpm_display_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function bpm_fft_display_Callback(hObject, eventdata, handles)

function bpm_fft_display_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function anomaly_display_Callback(hObject, eventdata, handles)

function anomaly_display_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ectopic_display_Callback(hObject, eventdata, handles)

function ectopic_display_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
