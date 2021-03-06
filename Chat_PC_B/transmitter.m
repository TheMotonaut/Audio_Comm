% COMPLETE THE TRANSMITTER!

% pack = message to be transmitted (consists of 432 bits from the GUI, always!)
% fc = carrier frequency



function transmitter(pack, fc)

constellation = [-1-i, -1+i, 1+i, 1-i]/sqrt(2);
M = length(constellation);                                   
bpsymb = log2(M);    

fs = 44000; %sampling frequency
Tsamp = 1/fs;
fc = 4000;
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

m = buffer(pack, bpsymb)';                         % Group bits into bits per symbol
m_idx = bi2de(m, 'left-msb')'+1;  
x = constellation(m_idx);


x_upsample = upsample(x, fsfd);

s = conv(y, x_upsample);

bq = [+1 +1 +1 +1 +1 -1 -1 +1 +1 -1 +1 -1 +1];
bq = resample(bq, 100, 1);
bq = bq.*cos(fc.*Tsamp.*(1:length(bq)));

I = sqrt(2).*real(s).*cos(2.*pi.*fc.*Tsamp.*(0:length(s)-1));
Q = sqrt(2).*imag(s).*sin(2.*pi.*fc.*Tsamp.*(0:length(s)-1));
tx = I + Q;
tx2 = tx;
tx = tx/max(abs(tx));
tx = [bq tx];
N0 = 70;
rx = awgn(tx, N0);
rx2 = awgn(tx2, N0);
correl = xcorr(rx, bq)
[M, I] = max(correl)

trimmed = rx(I-length(correl)/2+13:length(rx)-1)
%}

%figure(17)
%plot(tx)
%scatterplot(tx)
%sound(tx, freqs);
player = audioplayer(tx, fs);       %create an audioplayer object to play the noise at a given sampling frequency
playblocking(player); % Play the noise 

