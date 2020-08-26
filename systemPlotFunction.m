clear;
clc;
plot(-1.1861, 0 , 'kx','MarkerSize', 20);
hold on;
plot(-0.6469, 0.7277 , 'kx','MarkerSize', 20);
hold on;
plot(-0.6469, 0.7277 , 'kx','MarkerSize', 20);
hold on;
plot(-0.1613, 0.7701 , 'kx','MarkerSize', 20);
hold on;
plot(-0.1613, -0.7701 , 'kx','MarkerSize', 20);
hold on;
plot(0, 1.0515 , 'ko','MarkerSize', 20);
hold on;
plot(0, -1.0515 , 'ko','MarkerSize', 20);
hold on;
plot(0, 1.7013 , 'ko','MarkerSize', 20);
hold on;
plot(0, -1.7013 , 'ko','MarkerSize', 20);
grid on;
title('System Zeros and Poles');
axis([-1.5 0.5 -2 2]);