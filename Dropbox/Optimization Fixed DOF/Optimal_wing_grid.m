
%%% Subroutine to plot wing edges %%%%%%%

clear all
close all
clc


% linear
x1=load('Gridx_P3.dat');
y1=load('Gridy_P3.dat');
z1=load('Gridz_P3.dat');

[nrp1,ncp1] = size(x1);

x2=load('Gridx_P3-1.dat');
y2=load('Gridy_P3-1.dat');
z2=load('Gridz_P3-1.dat');

[nrp2,ncp2] = size(x2);

% % quadratic
% x2=load('Gridx_P2.dat');
% y2=load('Gridy_P2.dat');
% z2=load('Gridz_P2.dat');
% 
% [nrp2,ncp2] = size(x2);
% 
% % cubic
% x3=load('Gridx_P3.dat');
% y3=load('Gridy_P3.dat');
% z3=load('Gridz_P3.dat');
% 
% [nrp3,ncp3] = size(x3);
% 
% % quartic
% x4=load('Gridx_P4.dat');
% y4=load('Gridy_P4.dat');
% z4=load('Gridz_P4.dat');

% [nrp4,ncp4] = size(x4);

for i = 1:nrp1
plot3 ( x1(i, 1:ncp1), y1(i, 1:ncp1), z1(i, 1:ncp1), 'r', 'linewidth', 2 )
hold on
end

for j = 1:ncp1
    plot3 ( x1(1:nrp1, j), y1(1:nrp1, j), z1(1:nrp1, j), 'r', 'linewidth', 2 )
end

hold on

for i = 1:nrp2
plot3 ( x2(i, 1:ncp2), y2(i, 1:ncp2), z2(i, 1:ncp2), ':b', 'linewidth', 2 )
hold on
end

for j = 1:ncp2
    plot3 ( x2(1:nrp2, j), y2(1:nrp2, j), z2(1:nrp2, j), ':b', 'linewidth', 2 )
end

axis equal