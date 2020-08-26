% Author: Vasilis Mpletsos
% AEM: 8687
% FFT Low Pass Inverse Chebyshev
% 24/08/2020

close all;
clear all;
clc

%% Initiallization of values
target = 5; %dB target gain
aem = [8 6 8 7];
m = 2; 
f_p = (4+m)*(10^(3));
f_s = f_p/2.6;
a_min = 24 + aem(3)*(6/9);
a_max = 0.5 + a_min/36;

%% Calculate Filter
w_p = 2*pi*f_p;
w_s = 2*pi*f_s;
W_s = w_p/w_s;
W_p = 1;
n  = log10((10^(a_min/10)-1)/(10^(a_max/10)-1))/(2*log10(W_s));
n = ceil(n);
w_hp = 1/((10^(a_max/10)-1)^(1/(2*n)));
w0 = w_p/ w_hp;
s_1 = -1;
s_2 = -0.809 + 0.587 * 1i; 
w_2 = -0.809 - 0.587 * 1i;
s_3 = -0.309 + 0.951 * 1i;
w_3 = -0.309 - 0.951 * 1i;
Q1 = 0.5;
Q2 = 0.618;
Q3 = 1.618;

%% Unit 1
R1 = 1000;
km = 1000;
kf = w0;
C1 = 1/(km*kf);

%% Unit 2
R1_2 = 2*Q2;
R2_2 = 1/(2*Q2);
C1_2 = 1;
C2_2 = 1;
k = 1;
km = C2_2*(10^(8)/kf);
R1_2 = km * R1_2;
R2_2 = km * R2_2;
C1_2 = C1_2*(1/(km*kf));
C2_2 = C2_2*(1/(km*kf));

%% Unit 3
R1_3 = 2*Q3;
R2_3 = 1/(2*Q3);
C1_3 = 1;
C2_3 = 1;
k = 1;
km = C2_3*(10^(8)/kf);
R1_3 = km * R1_3;
R2_3 = km * R2_3;
C1_3 = C1_3*(1/(km*kf));
C2_3 = C2_3*(1/(km*kf));

%% Transfer Functions

Gain= 10^(target/20);
a = Gain/k;
R1 = 1000;
R2 = R1*a;

T1 = tf([1 0],[1 ( w0 )]);
T2 = tf([1 0 0],[1 ( w0 / Q2 ) (w0)^2]);
T3 = tf([1 0 0],[1 ( w0 / Q3 ) (w0)^2]);

T = T1*T2*T3;
T = T*a;
Tinv = inv(T);

f_gain = 40000;

plot_transfer_function(T1,[f_p f_s]);
plot_transfer_function(T2,[f_p f_s]);
plot_transfer_function(T3,[f_p f_s]);
plot_transfer_function(T,[f_p f_s f_gain]);
plot_transfer_function(Tinv,[f_p f_s f_gain]);

%% Calculate Added Signal Frequencies
w_1 = 0.5*w_s;
f_1 = w_1/(2*pi);
w_2 = 0.8*w_s;
f_2 = w_2/(2*pi);
w_3 = 1.2*w_p;
f_3 = w_3/(2*pi);
w_4 = 2.4*w_p;
f_4 = w_4/(2*pi);
w_5 = 3.5*w_p;
f_5 = w_5/(2*pi);

%% FFT Section
FFT_HP_B_Vmpletsos;
