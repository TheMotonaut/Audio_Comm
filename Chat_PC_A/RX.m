clear all; clc

freqs = 44000;
Tsamp = 1/freqs;
tau = 1/440; % 
alpha = 0.35;
span = 6;
BW = (1+alpha)/(2*tau);

y = rtrcpuls(alpha, tau, freqs, span);

constellation = [-1-i, -1+i, 1+i, 1-i]/sqrt(2);
%constellation = [(3 + 3i), (3 + 1i), (1 + 3i), (1 + 1i), (-3 + 3i), (-3 + 1i), (-1 + 3i), (-1 + 1i), (-3 - 3i), (-3 - 1i), (-1 - 3i), (-1 - 1i), (3 - 3i), (3 - 1i), (1 - 3i), (1 - 1i)];

deviceReader = audioDeviceReader(sampleRateValue)
deviceReader();
rx_temp = rx;
%rx = tx;
%plot(real(s))
%title('16QAM Real part of S before fc: 10dB SNR')


I_rx = sqrt(2).*rx.*cos(2.*pi.*fc.*Tsamp.*(0:length(s)-1));
Q_rx = 1i.*sqrt(2).*rx.*sin(2.*pi.*fc.*Tsamp.*(0:length(s)-1));
%subplot(2,2,3)
rx = lowpass(I_rx + Q_rx, fc/2, freqs);
%rx = I_rx + Q_rx;
%plot(real(rx))
%title('16QAM Real part of rx after fc and lowpass filter: 10dB SNR')

MF = fliplr(conj(y));

MF_out = conv(MF, rx);

%subplot(2,2,4)
%plot(real(MF_out))
%title('16QAM RX after MF: 10dB SNR')
MF_out = MF_out(2*span*freqs*tau  :  end-2*span*freqs*tau);

rx_vec = MF_out(1:freqs*tau:end);

L = length(s);
Y = fft(rx_vec);

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = freqs*(0:(L/2))/L;
%figure(3)
plot(f, P1)
%scatterplot(rx_vec);
data_out = {};
for point = rx_vec
    temp = 10000;
    temp2 = 0;
    i = 1;
    for const = constellation
        if(abs(point - const) < temp)
            temp = abs(point - const);
            temp2 = i;
        end
        i = i + 1;
    end
    
    data_out =[data_out, temp2];
end
data_out = cell2mat(data_out);
  
if all(data == data_out)
    disp('ok!')
else
    disp('not ok')
end

%sound(rx_temp, freqs)
%eyediagram(tx, freqs*tau)
%plot(real(s)))
%scatterplot(s)
%yediagram(s, freqs*tau)