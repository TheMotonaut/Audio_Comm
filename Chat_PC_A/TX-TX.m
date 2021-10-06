clear all; clc

freqs = 48000;
Tsamp = 1/freqs;
tau = 1/480; % 
alpha = 0.35;
span = 6;
BW = (1+alpha)/(2*tau);

y = rtrcpuls(alpha, tau, freqs, span);

constellation = [-1-i, -1+i, 1+i, 1-i]/sqrt(2);
constellation = [(3 + 3i), (3 + 1i), (1 + 3i), (1 + 1i), (-3 + 3i), (-3 + 1i), (-1 + 3i), (-1 + 1i), (-3 - 3i), (-3 - 1i), (-1 - 3i), (-1 - 1i), (3 - 3i), (3 - 1i), (1 - 3i), (1 - 1i)];

N = 10000;
data = randsrc(1,N,[0 1 2 3])+1;1
data = randsrc(1,N,[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15])+1;1

x = constellation(data);
x_upsample = upsample(x, freqs*tau);

s = conv(y, x_upsample);
fc = 5000;
IQ = s.*exp(-1i*2*pi*fc*(0:length(s)-1)*Tsamp);

I = sqrt(2).*real(s).*cos(2.*pi.*fc.*Tsamp.*(0:length(s)-1));
Q = sqrt(2).*imag(s).*sin(2.*pi.*fc.*Tsamp.*(0:length(s)-1));

tx = I + Q;
tx = tx/max(abs(tx));

sound(tx, freqs);
%figure(1)
%subplot(2,2,1)
%plot(tx)
%title('16QAM Transmitted signal with IQ and fc: 10dB SNR')
