% COMPLETE THE TRANSMITTER!

% pack = message to be transmitted (consists of 432 bits from the GUI, always!)
% fc = carrier frequency



function transmitter(pack, fc)

constellation = [-1-i, -1+i, 1+i, 1-i]/sqrt(2);
M = length(constellation);                                   
bpsymb = log2(M);    

fs = 44000; %sampling frequency
Tsamp = 1/fs;
fc = 5000;
Rb = 440;
fsymb = Rb/bpsymb;
Tsymb = 1/fsymb;
fsfd = fs/fsymb;
alpha = 0.35;
span = 6;
BW = (1+alpha)/(2*Tsymb);



N = 432;

%pack = randsrc(1,N,[0 1]);  % or "randi(2,1,N)-1"  %
%pack = [1,1,1,1,1,1,1,1,0,0,1,1,0,0,1,1,1,1,0,0,0,0,1,0,1,0,1,0,1,0,1,0,1,0,1,1,0,1,1,0,0,1,0,0,1,1,0,0,0,1,1,1,0,1,0,1,1,0,1,0,1,1,0,0,0,1,1,1,1,1,0,1,0,1,0,0,1,0,0,1,0,0,0,1,0,1,1,0,0,0,0,0,0,0,1,0,1,1,0,0];
y = rtrcpuls(alpha, Tsymb, fs, span);   


                              % Number of symbols in the constellation

m = buffer(pack, bpsymb)';                         % Group bits into bits per symbol
m_idx = bi2de(m, 'left-msb')'+1;  


test = pack;
save('file.mat', 'test')

%disp(pack)
%disp(size(pack))
%disp(size(m_idx))

x = constellation(m_idx);

%disp(length(x))
x_upsample = upsample(x, fsfd);

s = conv(y, x_upsample);

%scatterplot(awgn(x, 20))
figure(13)
bq = [+1 +1 +1 +1 +1 -1 -1 +1 +1 -1 +1 -1 +1];
%subplot(1,2,1)
%plot(bq)
bq = resample(bq, 100, 1);
%subplot(1,2,2)
plot(bq)


bq = bq.*cos(fc.*Tsamp.*(1:length(bq)));
%plot(bq)

I = sqrt(2).*real(s).*cos(2.*pi.*fc.*Tsamp.*(0:length(s)-1));
Q = sqrt(2).*imag(s).*sin(2.*pi.*fc.*Tsamp.*(0:length(s)-1));

tx = I + Q;
tx2 = tx;
tx = tx/max(abs(tx));
tx = [bq tx];

%figure(26)
%plot(tx)

%disp(pack)

%figure(19)
%plot(real(tx))

%tx = [zeros(1,500) tx];
%{

figure(20)
%plot(real(tx))
%plot(xcorr(tx, bq))

N0 = 70;
rx = awgn(tx, N0);
rx2 = awgn(tx2, N0);

subplot(2,2,1)
plot(xcorr(rx, bq))
subplot(2,2,2)
plot(rx)

correl = xcorr(rx, bq)
[M, I] = max(correl)

trimmed = rx(I-length(correl)/2+13:length(rx)-1)
subplot(2,2,3)
plot(trimmed)

subplot(2,2,4)
plot(rx2)
%}

%figure(17)
%plot(tx)
%scatterplot(tx)
%sound(tx, freqs);
%tx_signal = randn(1,1e4); %vector of noise
player = audioplayer(tx, fs);       %create an audioplayer object to play the noise at a given sampling frequency
playblocking(player); % Play the noise 
%subplot(1,2,1)
%plot(real(x))
%subplot(1,2,2)
%plot(real(x_upsample))
