
clear all
close all
clc
rho = 1.225;
V_ref = 10;

LBL = load('LOADS_BL.dat');
LOPT = load('AEROLOADSWR_P3.dat');
%LOPT2 = load('LOADS_BL.dat');

phi = 45*cos(2*pi*LBL(:,1));

% Thrust
figure(1)
plot(phi(10:end),-LBL(10:end,2),'LineWidth',2)
hold on
plot(phi(10:end),-LOPT(10:end,2),'r','LineWidth',2)
%hold on
%plot(phi(10:end),-LOPT2(10:end,2),':k')
xlabel('\phi ','fontsize',14)
ylabel('T^*','fontsize',14)   
legend('Baseline case','Optimal case','fontsize',14)

% Lift
figure(2)
plot(phi(10:end),LBL(10:end,3),'LineWidth',2)
hold on
plot(phi(10:end),LOPT(10:end,3),'r','LineWidth',2)
%hold on
%plot(phi(10:end),LOPT2(10:end,3),':k')
xlabel('\phi ','fontsize',14)
ylabel('L^*','fontsize',14)   
legend('Baseline case','Optimal case','fontsize',14)

% Aerodynamic power
figure(3)
plot(phi(10:end),LBL(10:end,4)/(rho*V_ref^3),'LineWidth',2)
hold on
plot(phi(10:end),LOPT(10:end,4)/(rho*V_ref^3),'r','LineWidth',2)
%hold on
%plot(phi(10:end),LOPT2(10:end,4),':k')
xlabel('\phi ','fontsize',14)
ylabel('Aerodynamic power','fontsize',14)   
legend('Baseline case','Optimal case','fontsize',14)