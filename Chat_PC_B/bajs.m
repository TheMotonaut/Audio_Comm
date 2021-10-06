clear all; clc

constellation = [-1-i, -1+i, 1+i, 1-i];

fs = 44000;
fc = 5000;
Tsamp = 1/fs;
M = length(constellation);
bpsymb = log2(M);
Rb = 440;
fsymb = Rb/bpsymb;
Tsamp = 1/fs;
alpha = 0.35;
span = 6;
Tsymb = 1/fsymb;
BW = (1+alpha)/(2*Tsymb);
fsfd = fs/fsymb;

y = rtrcpuls(alpha, Tsymb, fs, span);


%constellation = [(3 + 3i), (3 + 1i), (1 + 3i), (1 + 1i), (-3 + 3i), (-3 + 1i), (-1 + 3i), (-1 + 1i), (-3 - 3i), (-3 - 1i), (-1 - 3i), (-1 - 1i), (3 - 3i), (3 - 1i), (1 - 3i), (1 - 1i)];

N = 216;
data = randsrc(1,N,[0 1 2 3])+1;
%data = randsrc(1,N,[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15])+1;

%bq = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1];
%bq = resample(bq, 100, 1);

%bq = bq.*cos(fc.*Tsamp.*(1:length(bq)));
%bq = bq/max(abs(bq));

load('bq.mat', '-mat');

figure(10)
plot(bq)
x = constellation(data);
x_upsample = upsample(x, fsfd);

s = conv(y, x_upsample);
fc = 5000;

Id = real(s).*cos(2.*pi.*fc.*Tsamp.*(0:length(s)-1));
Qd = imag(s).*sin(2.*pi.*fc.*Tsamp.*(0:length(s)-1));

tx = Id + Qd;
tx = tx/max(abs(tx));
txm = tx;
%tx = [zeros(1,5000) tx];
tx = [bq tx];

figure(11)
subplot(1,2,1)
plot(tx)
subplot(1,2,2)
plot(xcorr(tx, bq))

tx0 = [zeros(1,5000) tx];
tx0 = [tx0 zeros(1,5000)];
N0 = 5;
rx = awgn(tx0, N0);

rx = rx/max(abs(rx));

rx_temp = rx;
z = 6;
figure(4)
subplot(2,2,1)
plot(rx)
correl = xcorr(rx, bq);
[~, I] = max(correl);

rx = rx(I-length(correl)/2+1300:I-length(correl)/2+45560+1300);


subplot(2,2,2)
plot(correl)
subplot(2,2,3)
plot(rx)

I_rx = rx.*cos(2.*pi.*fc.*Tsamp.*(0:length(rx)-1));
Q_rx = 1i.*rx.*sin(2.*pi.*fc.*Tsamp.*(0:length(rx)-1));

rx = lowpass(I_rx + Q_rx, fc/5, fs);

MF = fliplr(conj(y));
MF_out = conv(MF, rx);
MF_out = MF_out(2*span*fs*Tsymb  :  end-2*span*fs*Tsymb);

eyediagram(MF_out, 200)

rx_vec = MF_out(1:fs*Tsymb:end);
rx_vec = rx_vec/max(abs(rx_vec));


scatterplot(rx_vec)

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
