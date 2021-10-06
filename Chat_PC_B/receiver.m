% RECEIVER 
%
% This is the receiver structure that you will have to complete.
% The function: receiver(fc) is a setup function for the receiver. Here,
% the audiorecorder object is initialized (see help audiorecorder or
% MATLAB's homepage for more information about the object).
% 
% The callback function audioTimerFcn() is a callback function that is
% triggered on a specified time interval (here it is determined in the
% setup function, by the variable time_value)
% 
% Your task is to extend this code to work in the project!
%%

function [audio_recorder] = receiver(fc)
fs = 44000; %sampling frequency
audio_recorder = audiorecorder(fs,24,1);% create the recorder

%attach callback function
time_value = 6; % how often the function should be called in seconds
set(audio_recorder,'TimerPeriod',time_value,'TimerFcn',@audioTimerFcn); % attach a function that should be called every second, the function that is called is specified below.

%ADD USER DATA FOR CALLBACK FUNCTION (DO NOT CHANGE THE NAMES OF THESE VARIABLES!)
audio_recorder.UserData.receive_complete = 0; % this is a flag that the while loop in the GUI will check
audio_recorder.UserData.pack  = []; %allocate for data package
audio_recorder.UserData.pwr_spect = []; %allocate for PSD
audio_recorder.UserData.const = []; %allocate for constellation
audio_recorder.UserData.eyed  = []; %allocate for eye diagram


record(audio_recorder); %start recording
end


% CALLBACK FUNCTION
% This function will be called every [time_value] seconds, where time_value
% is specified above. Note that, as stated in the project MEMO, all the
% fields: pwr_spect, eyed, const and pack need to be assigned if you want
% to get outputs in the GUI.

% So, let's see an example of where we have a pulse train as in Computer
% exercise 2 and let the GUI plot it. Note that we will just create 432
% random bits and hence, the GUI will not be able to decode the message but
% only to display the figures.
% Parameters in the example:
% f_s = 22050 [samples / second]
% R_s = 350 [symbols / second]
% fsfd = f_s/R_s [samples / symbol]
% a = 432 bits
% M = 4 (using QPSK as in computer exercise)

function audioTimerFcn(recObj, event, ~)

%-----------------------------------------------------------
% THE CODE BELOW IS BASED ON COMPUTER EX 5 AND EX 6:
%-----------------------------------------------------------
disp('Callback triggered')                                         % Number of samples per symbol (choose fs such that fsfd is an integer) [samples/symb]                          % create sinc pulse with span = 6

constellation = [-1-1i, -1+1i, 1+1i, 1-1i];
N = 432;
fs = 44000;
fc = 6000;
Tsamp = 1/fs;
M = length(constellation);
bpsymb = log2(M);
Rb = 440;
fsymb = Rb/bpsymb;
Tsamp = 1/fs;
tau = 1/440; % 
alpha = 0.35;
span = 6;
Tsymb = 1/fsymb;
BW = (1+alpha)/(2*tau);
fsfd = fs/fsymb;
%preamble = [1 1 1 -1 -1 1 -1];     % Length-7 Barker code
%disp('2')

bq = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1];
%cos(pi)
%disp('3')
bq = resample(bq, 100, 1);
bq = bq.*cos(fc.*Tsamp.*(1:length(bq)));

figure(40)
plot(bq)
%('1')
y = rtrcpuls(alpha, Tsymb, fs, span);

rx = getaudiodata(recObj).';

%rx = rx/max(abs(rx));

rx_temp = rx;
figure(39)
%subplot(2,2,1)
%plot(rx)
correl = xcorr(rx, bq);
%subplot(2,2,2)
%plot(correl)    
[~, I] = max(correl);
%disp(I)
%disp(length(rx))
rx = rx(I-length(correl)/2+1300:45560+I-length(correl)/2+1300);

%subplot(2,2,3)
%plot(rx)
I_rx = rx.*cos(2.*pi.*fc.*Tsamp.*(0:length(rx)-1));
Q_rx = 1i.*rx.*sin(2.*pi.*fc.*Tsamp.*(0:length(rx)-1));

rx = lowpass(I_rx + Q_rx, fc/10, fs);

MF = fliplr(conj(y));
MF_out = conv(MF, rx);
MF_out = MF_out(2*span*fs*Tsymb  :  end-2*span*fs*Tsymb);
rx_vec = MF_out(1:fs*Tsymb:end);
rx_vec = rx_vec/median(abs(rx_vec));
%subplot(2,2,4)
%plot(real(rx))
%disp(rx_vec)
% data_out = {};
% for point = rx_vec
%     temp = 10000;
%     temp2 = 0;
%     i = 1;
%     for const = constellation
%         if(abs(point - const) < temp)
%             temp = abs(point - const);
%             temp2 = i;
%         end
%         i = i + 1;
%         
%     end
%     data_out = [data_out, temp2];
% end

rx_vec = transpose(rx_vec);
constellation = transpose(constellation);

const_rep = repmat(constellation, [1, length(rx_vec)]);

[~, index] = min(transpose(abs(rx_vec - transpose(const_rep))));

test = [0;1;1;0;0;1;0;0;0;1;1;0;0;0;0;1;0;1;1;1;0;0;1;1;0;1;1;0;1;0;1;1;0;1;1;0;0;1;0;0;0;1;1;0;1;1;1;1;0;1;1;0;0;0;0;1;0;1;1;1;0;0;1;1;0;1;1;0;1;0;1;0;0;1;1;0;0;1;0;0;0;1;1;0;1;0;1;0;0;1;1;0;0;0;0;1;0;1;1;0;1;0;0;1;0;1;1;1;0;0;1;1;0;1;1;0;0;1;0;0;0;1;1;0;1;0;1;0;0;1;1;1;0;0;1;1;0;1;1;0;0;0;0;1;0;1;1;0;1;1;0;1;0;1;1;0;0;1;0;0;0;1;1;0;1;0;1;1;0;1;1;0;0;0;0;1;0;1;1;1;0;0;1;1;0;1;1;0;1;0;1;0;0;1;1;0;0;1;0;0;0;1;1;0;1;0;1;1;0;1;1;0;1;0;0;1;0;1;1;0;0;0;0;1;0;1;1;0;1;0;1;0;0;1;1;0;0;1;0;0;0;1;1;0;1;0;0;1;0;1;1;0;1;0;1;0;0;1;1;0;0;0;0;1;0;1;1;0;0;1;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;1;1;1;1;1;0;1;1;1;0;0;1;0;0;1;0;1;1;0;1;1;1;1;1;1;0;1;1;1;0;0];
%load('file.mat')
%data_out = cell2mat(data_out)-1;

%disp(sum(data_out))


pulse_train = MF_out;
x = rx_vec;
%data_out = de2bi(data_out, 'left-msb');
m = de2bi(index - 1, bpsymb, 'left-msb'); 
data_out = reshape(m', N, 1);
% disp(size(test))
% disp(size(data_out))
% data_out = data_out(:).';
% disp(data_out)
% figure(592)
% subplot(2,1,1)
% plot(test)
% subplot(2,1,2)
% plot(data_out)
if all(test == data_out.')
    disp('ok!')
else
    disp('not ok')
end
%disp(length(data_out))

%------------------------------------------------------------------------------
% HOW TO SAVE DATA FOR THE GUI
%   NOTE THAT THE EXAMPLE HERE IS ONLY USED TO SHOW HOW TO OUTPUT DATA
%------------------------------------------------------------------------------

% Step 1: save the estimated bits
recObj.UserData.pack = data_out;

% Step 2: save the sampled symbols
recObj.UserData.const = x;

% Step 3: provide the matched filter output for the eye diagram
recObj.UserData.eyed.r = pulse_train;
recObj.UserData.eyed.fsfd = fsfd;

% Step 4: Compute the PSD and save it. 
% !!!! NOTE !!!! the PSD should be computed on the BASE BAND signal BEFORE matched filtering
[pxx, f] = pwelch(pulse_train,1024,768,1024, fs); % note that pwr_spect.f will be normalized frequencies
f = fftshift(f); %shift to be centered around fs
f(1:length(f)/2) = f(1:length(f)/2) - fs; % center to be around zero
p = fftshift(10*log10(pxx/max(pxx))); % shift, normalize and convert PSD to dB
recObj.UserData.pwr_spect.f = f;
recObj.UserData.pwr_spect.p = p;

% In order to make the GUI look at the data, we need to set the
% receive_complete flag equal to 1:
recObj.UserData.receive_complete = 1; 
    
end
