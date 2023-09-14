clear 
clc; 
close all;
%% initializations
nUsers =2;
nBS = 3;
nTx = 2;
P_sp = 10*10^-2; %check 10 W/Gnats/s, for BW=10^6 
C_b = 30 ; %Mbps
P_bc = 10^(30/10)/1000; %[mW] BS power constraint

load ../data/testChannel.mat;
[cbv,cup,A] = BRB_main(channel,scale,nUsers,nBS,nTx,P_bc,P_sp,C_b) ;
plot(A(:,1))
hold on
plot(A(:,2));
plot(A(:,3));
legend("Upper bound", "Current best objective","Lower bound")
saveas(gcf, '../../results/ConvergencePlot.png')