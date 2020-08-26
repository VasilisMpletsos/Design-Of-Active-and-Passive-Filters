% Author: Vasilis Mpletsos
% AEM: 8687
% Band Pass Inverse Chebyshev
% 17/08/2020

close all
clear all
clc

%% Initiallization
target = 5; %dB target gain
aem = [8 6 8 7];
f0 = 1000;
f1 = 650 + 25*aem(4);
f2 = (f0^2)/f1;
D = 2.1 * (((f0^2)-f1^2)/f1);
f3 = (-D + sqrt((D^2)+(4*(f0^2))))/2;
f4 = (f0^2)/f3;
a_min = 35 - aem(3);
a_max = 0.4 + (aem(4)/36);
w0 = 2*pi*f0;
w1 = 2*pi*f1;
w2 = 2*pi*f2;
w3 = 2*pi*f3;
w4 = 2*pi*f4;

%% Filter Design
bw = w2-w1;
w_p = 1;
w_s = (w4-w3)/(w2-w1);
w_p_norm = w_p/w_s; 
w_s_norm = 1;
n = acosh(sqrt((((10^(a_min/10)) -1))/ ((10^(a_max/10))- 1)))/acosh(1/w_p_norm);
n = ceil(n);
e = (10^(a_min/10)-1)^(-1/2);
a = (1/n)*(asinh(1/e));
W_hp = 1/(cosh((1/n)*(acosh(1/e))));

%% Chebyshev for angles (+- 22.5 and +-67.5)
s_0 = -sinh(a)*cosd(22.5);
s_1 = -sinh(a)*cosd(67.5);
W_0 = cosh(a)*sind(22.5);
W_1 = cosh(a)*sind(67.5);
w_0 = sqrt( s_0^2 + W_0^2);
w_1 = sqrt( s_1^2 + W_1^2);
Q_0 = 1/(2*cos(atan(W_0/s_0)));
Q_1 = 1/(2*cos(atan(W_1/s_1)));
p1 = -sinh(a)*cosd(22.5)+(1i)*cosh(a)*sind(22.5);
p2 = -sinh(a)*cosd(-22.5)+(1i)*cosh(a)*sind(-22.5);
p3 = -sinh(a)*cosd(67.5)+(cosh(a)*sind(67.5))*1i;
p4 = -sinh(a)*cosd(-67.5)+(cosh(a)*sind(-67.5))*1i;

%% Inverse Chevyshev
winv_0 = (1/w_0)*(1/w_p_norm);
winv_1 = (1/w_1)*(1/w_p_norm);
S_0 = (winv_0/2)*(1/Q_0);
S_1 = (winv_1/2)*(1/Q_1);
Winv_0 = sqrt((winv_0^2)-(S_0^2));
Winv_1 = sqrt((winv_1^2)-(S_1^2));
wz_1 = sec(pi/(2*n));
wz_2 = sec((3*pi)/(2*n));
wz_1 = wz_1*w_s;
wz_2 = wz_2*w_s;

%% Geffe Section
qc = w0/bw;

% First 
C1 = S_0^2 + Winv_0^2;
D1 = (2*S_0)/qc;
E1 = 4 + (C1/(qc^2));
G1 = sqrt((E1^2)-(4*(D1^2)));
Q1 = (1/D1)*sqrt((E1+G1)/2);
k1 = (S_0*Q1)/qc;
W1 = k1 + sqrt((k1^2)-1);
w01 = (1/W1)*w0;
w02 = W1*w0;

% Second
C2 = S_1^2 + Winv_1^2;
D2 = (2*S_1)/qc;
E2 = 4 + (C2/(qc^2));
G2 = sqrt((E2^2)-(4*(D2^2)));
Q2 = (1/D2)*sqrt((E2+G2)/2);
k2 = (S_1*Q2)/qc;
W2 = k2 + sqrt((k2^2)-1);
w03 = (1/W2)*w0;
w04 = W2*w0;

% Transform imaginary zero wz_1
kimg_1 = 2 + ((wz_1^2)/(qc^2));
x1 = (kimg_1+sqrt((kimg_1^2)-4))/2;
wimgz_1 = sqrt(x1)*w0;
wimgz_2 = w0/sqrt(x1);

% Transform imaginary zero wz_2
kimg_2 = 2 + ((wz_2^2)/(qc^2));
x2 = (kimg_2+sqrt((kimg_2^2)-4))/2;
wimgz_3 = sqrt(x2)*w0;
wimgz_4 = w0/sqrt(x2);

%% Unit 1
w0 = w01;
wz = wimgz_1;
wz_norm = wz/w0;
w0 = 1;
fprintf('k1 should be larger than %d and smaller than 1 \n',w0^2/wz_norm^2);
k1 = 0.5; %I desided randomly k1 inside the legal margin
R1_1 = (2/((k1*(wz_norm^2))-1));
R2_1 = 1/(1-k1);
R3_1 = 0.5*((k1/(Q1^2))+k1*(wz_norm^2)-1);
R4_1 = 1/k1;
R5_1 = 1;
R6_1 = 1;
C1_1 = k1/(2*Q1);
C2_1 = 2*Q1;
k_1 = R5_1/(R5_1+R3_1);
kf = w01;
km = C1_1*(10^(7)/kf);
R1_1 = km * R1_1;
R2_1 = km * R2_1;
R3_1 = km * R3_1;
R4_1 = km * R4_1;
R5_1 = km * R5_1;
R6_1 = km * R6_1;
C1_1 = C1_1*(1/(km*kf));
C2_1 = C2_1*(1/(km*kf));


%% Unit 2

% BoctorHighPass(wimgz_2,w02,Q1);
% A Boctor HPN cannot be defined. You should use a High Pass Notch (Fig. 7.21).
WZ_2 =  w02/wimgz_2;
k2_1 = (WZ_2^2)-1;
k2_2 = ((2+k2_1)*(Q1^2))/((2+k2_1)*(Q1^2)+1);
R1_2 = 1;
R2_2 = (Q1^2)*((k2_1+2)^2);
R3_2 = 1;
R4_2 = (Q1^2)*(k2_1+2);
C2 = 1/((Q1)*(k2_1+2));
C2_2 = k2_1 * C2;
k_2 = k2_2 * (WZ_2^2);
kf = w02;
km = C2*((10^(7))/kf);
R1_2 = km * R1_2;
R2_2 = km * R2_2;
R3_2 = km * R3_2;
R4_2 = km * R4_2;
C2 = C2/(km*kf);
C2_2 = C2_2/(km*kf);

%% Unit 3
w0 = w03;
wz = wimgz_3;
wz_norm = wz/w0;
w0 = 1;
fprintf('k1 should be larger than %d and smaller than 1 \n',w0^2/wz_norm^2);
k1 = 0.5; %I desided randomly k1 inside the legal margin
R1_3 = (2/((k1*(wz_norm^2))-1));
R2_3 = 1/(1-k1);
R3_3= 0.5*((k1/(Q2^2))+k1*(wz_norm^2)-1);
R4_3 = 1/k1;
R5_3 = 1;
R6_3 = 1;
C1_3 = k1/(2*Q2);
C2_3 = 2*Q2;
k_3 = R5_3/(R5_3+R3_3);
kf = w03;
km = C1_3*((10^(7)/kf));
R1_3 = km * R1_3;
R2_3 = km * R2_3;
R3_3 = km * R3_3;
R4_3 = km * R4_3;
R5_3 = km * R5_3;
R6_3 = km * R6_3;
C1_3 = C1_3*(1/(km*kf));
C2_3 = C2_3*(1/(km*kf));

%% Unit 4

% BoctorHighPass(wimgz_4,w04,Q2);
% A Boctor HPN cannot be defined. You should use a High Pass Notch (Fig. 7.21).
WZ_4 =  w04/wimgz_4;
k4_1 = (WZ_4^2)-1;
k4_2 = ((2+k4_1)*(Q2^2))/((2+k4_1)*(Q2^2)+1);
R1_4 = 1;
R2_4 = (Q2^2)*(k4_1+2)^2;
R3_4 = 1;
R4_4 = (Q2^2)*(k4_1+2);
C4 = 1/((Q2)*(k4_1+2));
C4_1 = k4_1 * C4;
k_4 = k4_2 * (WZ_4^2);
kf = w04;
km = C4*((10^(7))/kf);
R1_4 = km * R1_4;
R2_4 = km * R2_4;
R3_4 = km * R3_4;
R4_4 = km * R4_4;
C4 = C4/(km*kf);
C4_1 = C4_1/(km*kf);

%% Transfer Functions
w0 = 2*pi*f0;

T1 = tf([k_1*1, 0, k_1*(wimgz_1^2)], [1, w01/Q1, w01^2]);
T2 = tf([k_2*1, 0, k_2*(wimgz_2^2)], [1, w02/Q1, w02^2]);
T3 = tf([k_3*1, 0, k_3*(wimgz_3^2)], [1, w03/Q2, w03^2]);
T4 = tf([k_4*1, 0, k_4*(wimgz_4^2)], [1, w04/Q2, w04^2]);

T = T1*T2*T3*T4;
Gain = k_1*k_2*k_3*k_4;

Tbp_1 = k_1*((wimgz_1^2-w0^2)/(sqrt((w01^2-w0^2)^2 + ((w01*w0)/Q1)^2)));
Tbp_2 = k_2*((wimgz_2^2-w0^2)/(sqrt((w02^2-w0^2)^2 + ((w02*w0)/Q1)^2)));
Tbp_3 = k_3*((wimgz_3^2-w0^2)/(sqrt((w03^2-w0^2)^2 + ((w03*w0)/Q2)^2)));
Tbp_4 = k_4*((wimgz_4^2-w0^2)/(sqrt((w04^2-w0^2)^2 + ((w04*w0)/Q2)^2)));
Tbp = Tbp_1*Tbp_2*Tbp_3*Tbp_4;

a = (10^(target/20))/Tbp;
R1 = 10^4;
R2 = R1 * a;
T = a*T; 
Tinv = inv(T);

%% Plot Section
plot_transfer_function(T1, [f1 f2 f3 f4])
plot_transfer_function(T2, [f1 f2 f3 f4])
plot_transfer_function(T3, [f1 f2 f3 f4])
plot_transfer_function(T4, [f1 f2 f3 f4])
ltiview({'bodemag'}, T1, T2, T3, T4, T)
plot_transfer_function(T, [f1 f2 f3 f4])
plot_transfer_function(Tinv, [f1 f2 f3  f4])

%% Calculate Added Signal Frequencies
w_1 = w0 - ((w0-w1)/3);
f_1 = w_1/(2*pi);
w_2 = w0 + ((w0+w1)/4);
f_2 = w_2/(2*pi);
w_3 = 0.5*w3;
f_3 = w_3/(2*pi);
w_4 = 2.4*w4;
f_4 = w_4/(2*pi);
w_5 = 3*w4;
f_5 = w_5/(2*pi);

%% Final Step for FFT
FFT_BP_IC_Vmpletsos
