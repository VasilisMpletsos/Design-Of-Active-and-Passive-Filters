% Author: Vasilis Mpletsos
% AEM: 8687
% FFT Low Pass Inverse Chebyshev
% 12/08/2020

close all;
clear all;
clc

%% Initiallization of values
aem = [8 6 8 7];
m = 1; 
f_p = 1.1*(3+m)*(10^(3)); %*1000 to make it kilo hertz
f_s = 1.7*f_p;
a_min = 25 + (max(1,aem(3))-5)*(3/4);
a_max = 0.5 + (max(1,aem(4))-5)/16;

%% Calculate Filter
w_p = 2*pi*f_p;
w_s = 2*pi*f_s;

%% Test Section of Example at PDF 9
%w_p = 1000;
%w_s = 1400;
%a_max = 0.25;
%a_min = 18;
%f_p = w_p/(2*pi);
%f_s = w_s/(2*pi);
capacitor = 1; %This is so i can play between examples 0.1uF and mine 0.01uF
target = 5/20; %This is so i can toggle between example and my gain target 


%% Proper Continue
w_p_norm = w_p/w_s;
w_s_norm = w_s/w_s;
n = acosh(sqrt((((10^(a_min/10)) -1))/ ((10^(a_max/10))- 1)))/acosh(1/w_p_norm);
n = ceil(n);

%% Find half of Power
e = (10^(a_min/10)-1)^(-1/2);
W_hp = 1/(cosh((1/n)*(acosh(1/e))));
a = (1/n)*(asinh(1/e));

%% Chebyshev for angles (0 and +- 36 and =-72)
s_0 = -sinh(a)*cosd(0);
s_1 = -sinh(a)*cosd(36);
s_2 = -sinh(a)*cosd(72);
W_0 = cosh(a)*sind(0);
W_1 = cosh(a)*sind(36);
W_2 = cosh(a)*sind(72);
w_0 = sqrt( s_0^2 + W_0^2);
w_1 = sqrt( s_1^2 + W_1^2);
w_2 = sqrt( s_2^2 + W_2^2);
Q_0 = 1/(2*cos(atan(W_0/s_0)));
Q_1 = 1/(2*cos(atan(W_1/s_1)));
Q_2 = 1/(2*cos(atan(W_2/s_2)));

%% Inverse Chebyshev
winv_0 = 1/w_0;
winv_1 = 1/w_1;
winv_2 = 1/w_2;
w_z_0 = sec((1*pi)/(2*n));
w_z_1 = sec((3*pi)/(2*n));
w_z_2 = sec((5*pi)/(2*n));

%% Unit 1
W_Z = w_z_0/winv_2;
if W_Z>1
    fprintf('*** Passed W_Z on Unit 1\n');
end
C_1 = 1/(2*Q_2);
R2_1 = 4*(Q_2^2);
R1_1 = 1;
R4_1 = 1;
R5_1 = R2_1/((W_Z^2)-1);
R3_1 = (W_Z^2)/(2*(Q_2^2));
kH_1 = 1/(R3_1 + 1);
kL_1 = kH_1*(W_Z^2);
kf_1 = w_s*winv_2;
km_1 = C_1/(kf_1*(10^(-8))*capacitor); %We want 0.01uF
C_1 = 0.01*10^(-6)*capacitor;
R1_1 = R1_1*km_1;
R2_1 = R2_1*km_1;
R3_1 = R3_1*km_1;
R4_1 = R4_1*km_1;
R5_1 = R5_1*km_1;

%% Unit 2
W_Z = w_z_1/winv_1;
if W_Z>1
    fprintf('*** Passed W_Z on Unit 2\n');
end
C_2 = 1/(2*Q_1);
R2_2 = 4*Q_1^2;
R1_2 = 1;
R4_2 = 1;
R5_2 = R2_2/((W_Z^2)-1);
R3_2 = W_Z^2/(2*(Q_1^2));
kH_2 = 1/(R3_2 + 1);
kL_2 = kH_2*(W_Z^2);
kf_2 = w_s*winv_1;
km_2 = C_2/(kf_2*(10^(-8))*capacitor); %We want 0.01uF
C_2 = 0.01*10^(-6)*capacitor;
R1_2 = R1_2*km_2;
R2_2 = R2_2*km_2;
R3_2 = R3_2*km_2;
R4_2 = R4_2*km_2;
R5_2 = R5_2*km_2;

%% Unit 3
R_3 = 1/winv_0;
kf_3 = w_s;
km_3 = 1/(kf_3*(10^(-8))*capacitor); %We want 0.01uF
R_3 = R_3*km_3;
C_3 = 10^(-8)*capacitor;

%% Gain
gain = kL_1*kL_2;
gain_db = 20*log10(gain);
goal = 10^target; %Due to my AEM i have to get 5dB gain
if gain == goal 
   fprintf('*** Passed Gain Requirements'); 
   k_add = 1;
else
    fprintf('*** Failed: We want 5dB and we have %.3fdB\n',gain_db);
    k_add = goal/gain; %%We have to add this!
    R1 = 10000;
    R2 = k_add * R1;
end

%`% Transfer Functions
w_z_0_correct = w_z_0 * w_s;
w_z_1_correct = w_z_1 * w_s;
winv_0_correct = winv_0 * w_s;
winv_1_correct = winv_1 * w_s; 
winv_2_correct = winv_2 * w_s; 

T1 = tf( [kH_1 0 kH_1*((w_z_0_correct)^2)] , [1 ((winv_2_correct)/Q_2) ((winv_2_correct)^2)] );
T2 = tf( [kH_2 0 kH_2*((w_z_1_correct)^2)] , [1 ((winv_1_correct)/Q_1) ((winv_1_correct)^2)] );
T3 = tf( [(winv_0_correct)] , [1 winv_0_correct]);

T = T1*T2*T3;
T = T*k_add;
inv_T = inv(T);

plot_transfer_function( T1, [0 f_p f_s]);
plot_transfer_function( T2, [0 f_p f_s]);
plot_transfer_function( T3, [0 f_p f_s]);
plot_transfer_function( T, [0 f_p f_s]);
plot_transfer_function(inv_T, [0 f_p f_s]);
ltiview({'bodemag'},T1,T2,T3,T);

%% Final Step for FFT
FFT_LP_IC_Vmpletsos
