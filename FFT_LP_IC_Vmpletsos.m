% Author: Vasilis Mpletsos
% AEM: 8687
% FFT Low Pass Inverse Chebyshev
% 16/08/2020


%% Initiallization
T_time = 10*(1/2000); %Generate a 10 periods of a sawtooth wave 
Fs = 1000000; %Sampling frequency 1 MHz
dt = 1/Fs;
t = 0:dt:T_time-dt;
x = 0.5*(sawtooth(2*pi*2000*t)+1); %2000 sawtooth

%% Plot Input
figure;
plot(t,x);
axis([0 inf -0.5 1.5]);
title('Input signal');

%% Plot Output
y = lsim(T,x,t);
figure;
plot(t,y);
title('Output signal');

%% Calculate FFT on Input
Xf=fft(x);
L=length(x);
P2 = abs(Xf/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
figure;
plot(f,P1);
axis([0.01 20000 0 inf]);
title('Input signal');
xlabel('(Hz)');
ylabel('(Vin)');

%% Calculate FFT on Output
Yf=fft(y);
L=length(y);
P2 = abs(Yf/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
figure;
plot(f,P1);
axis([0.01 20000 0 inf]);
title('Output signal');
xlabel('(Hz)');
ylabel('(Vout)');