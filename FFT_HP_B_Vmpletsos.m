% Author: Vasilis Mpletsos
% AEM: 8687
% FFT Band Pass Inverse Chebyshev
% 25/08/2020

%% Initiallization
Fs = 10^6;
Ts = 1/Fs;
t_ol = 3*(10^(-2));
t = 0:Ts:t_ol-Ts;
input = cos(w_1*t)+0.6*cos(w_2*t)+cos(w_3*t)+0.8*cos(w_4*t)+0.4*cos(w_5*t);

%% Plot Input
figure;
plot(t,input);
axis([0 0.003 -5 5]);
title('Input Signal');

%% Plot Output
output = lsim(T,input,t);
figure;
plot(t,output);
axis([0 0.003 -5 5]);
title('Output Signal');

%% Plot Together
figure;
plot(t,input);
hold on;
plot(t,output);
hold off;
axis([0 0.003 -5 5]);
legend({'Input','Output'},'Location','southwest')
title('Input and Output Signal');


%% Calculate FFT on Input
FFTinput = fft(input);
L = length(input);
P2 = abs(FFTinput/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
figure;
plot(f,P1);
axis([0.01 30000 0 1.2]);
title('Input FFT');


%% Calculate FFT on Output
FFToutput = fft(output);
L = length(output);
P2 = abs(FFToutput/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
figure;
plot(f,P1);
axis([0.01 30000 0 2]);
title('Output FFT');
