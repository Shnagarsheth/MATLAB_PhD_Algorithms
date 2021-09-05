%% --------Tunnel-diode oscillator fitlering eqautions------------%%
% Created by: Shaival Nagarsheth
% Contains - SDEs, unpertubed ODEs, EKF, Second-order filtering equaitons
% and third-order filtering equations
% Version & Date: 06-08-2020 (updated)
%%---------------------------------------------------------------%%

%% ------------Initialization------------------
clear all;
clc;
N=10000;
T=0.35e-5;
dt=T/N;
R=0.5;
alpha=2.5;
L=1e-6;
C=1e-9;
V=0.4;
% R=1500;
% alpha=0.0001;
% L=0.000005;
% C=0.00000000002;
% V=1.2;
x1(1)=0.131;
x2(1)=0.00004;


ux1(1)=0.131;
ux2(1)=0.00004;

Fmx1(1)=0.131;
Fmx2(1)=0.00004;
FP11(1)=0;
FP22(1)=0;
FP12(1)=0;
FP21(1)=0;

Smx1(1)=0.131;
Smx2(1)=0.00004;
SP11(1)=0;
SP22(1)=0;
SP12(1)=0;
SP21(1)=0;

Tmx1(1)=0.131;
Tmx2(1)=0.00004;
TP11(1)=0;
TP22(1)=0;
TP12(1)=0;
TP21(1)=0;

%% ------------------Main Loop-------------------
for k=1:N
    dW1t=sqrt(dt)*randn;
    dW2t=sqrt(dt)*randn;
    dW3t=sqrt(dt)*randn;
    
    %% SDEs    
    x1(k+1)=x1(k)+(((1/C)*(-(x1(k))^3+1.5*(x1(k))^2-0.6*(x1(k))+x2(k)))*dt)+(alpha*dW3t);
    x2(k+1)=x2(k)+(((1/L)*(-x1(k)-R*x2(k)+V))*dt)+(alpha*(dW1t+dW2t));
    
    %% Unperturbed
    ux1(k+1)=ux1(k)+(((1/C)*(-(ux1(k))^3+1.5*(ux1(k))^2-0.6*(ux1(k))+ux2(k)))*dt);
    ux2(k+1)=ux2(k)+(((1/L)*(-ux1(k)-R*ux2(k)+V))*dt);
    
    %% EKF - mean
    Fmx1(k+1)=Fmx1(k)+((1/C)*(((-Fmx1(k))^3+1.5*(Fmx1(k))^2-0.6*(Fmx1(k)))+Fmx2(k))*dt);
    Fmx2(k+1)=Fmx2(k)+((1/L)*(-Fmx1(k)-R*Fmx2(k)+V)*dt);
    
    %% EKF-Variance    
    FP11(k+1)=FP11(k)+((((2/C)*(FP11(k)*(-3*(Fmx1(k))^2+3*(Fmx1(k))-0.6))+((2/C)*(FP12(k))))+(alpha^2))*dt);    
    FP12(k+1)=FP12(k)+((((1/L)*(-FP11(k)))+((1/L)*(-R*FP12(k)))+((1/C)*(FP12(k)*(-3*(Fmx1(k))^2+(3*(Fmx1(k)))-0.6)))+((1/C)*(FP22(k))))*dt);    
    FP21(k+1)=FP12(k+1);
    FP22(k+1)=FP22(k)+((((2/L)*((-FP12(k))+(-R*FP22(k))))+(2*alpha^2))*dt);
    
    %% Second-order Filter mean
    Smx1(k+1)=Smx1(k)+(((1/C)*((((-Smx1(k))^3+1.5*(Smx1(k))^2-0.6*(Smx1(k)))+Smx2(k))+(0.5*(SP11(k))*(-6*(Smx1(k))+(3)))))*dt);
    Smx2(k+1)=Smx2(k)+((1/L)*(-Smx1(k)-R*Smx2(k)+V)*dt);
    
    %% Second-order Filter variance
    SP11(k+1)=SP11(k)+((((2/C)*(SP11(k)*(-3*(Smx1(k))^2+3*(Smx1(k))-0.6))+((2/C)*(SP12(k))))+(alpha^2))*dt);
    SP12(k+1)=SP12(k)+((((1/L)*(-SP11(k)))+((1/L)*(-R*SP12(k)))+((1/C)*(SP12(k)*(-3*(Smx1(k))^2+(3*(Smx1(k)))-0.6)))+((1/C)*(SP22(k))))*dt);
    SP21(k+1)=SP12(k+1);
    SP22(k+1)=SP22(k)+((((2/L)*((-SP12(k))+(-R*SP22(k))))+(2*alpha^2))*dt);
    
    %% Third-order filer mean
    Tmx1(k+1)=Tmx1(k)+(((1/C)*((((-Tmx1(k))^3+1.5*(Tmx1(k))^2-0.6*(Tmx1(k)))+Tmx2(k))+(0.5*(TP11(k))*(-6*(Tmx1(k))+(3)))))*dt);
    Tmx2(k+1)=Tmx2(k)+((1/L)*(-Tmx1(k)-R*Tmx2(k)+V)*dt);    
    
    %% Third-order fitler variance
    TP11(k+1)=TP11(k)+((((2/C)*(TP11(k)*(-3*(Tmx1(k))^2+3*(Tmx1(k))-0.6))+((2/C)*(TP12(k)))+((-6/C)*(TP11(k)*TP11(k))))+(alpha^2))*dt);
    TP12(k+1)=TP12(k)+((((1/L)*(-TP11(k)))+((1/L)*(-R*TP12(k)))+((1/C)*(TP12(k)*(-3*(Tmx1(k))^2+(3*(Tmx1(k)))-0.6)))+((1/C)*(TP22(k)))+((-3/C)*(TP12(k)*TP11(k))))*dt);
    TP21(k+1)=TP12(k+1);
    TP22(k+1)=TP22(k)+((((2/L)*((-TP12(k))+(-R*TP22(k))))+(2*alpha^2))*dt);    
end

t=0:dt:T;

%% Plots for mean and variance

% plot(t,x1);
% figure;
% plot(t,x2);
% figure;
% plot(t,x1,'r',t,ux1,'b');
% figure;
% plot(t,x2,'r',t,ux2,'b');
% figure;
% plot(t,ux1,'r',t,Fmx1,'b');
% figure;
% plot(t,ux2,'r',t,Fmx2,'b');

%usefull
% figure;
% plot(x2,x1);
% figure;
% plot(ux2,ux1);

% figure;
% plot(t,Fmx1);
% figure;
% plot(t,Fmx2);
% figure;
% plot(t,Smx1);
% figure;
% plot(t,Smx2);
% figure;
% plot(t,FP11,'r');
% figure;
% plot(t,FP22,'g');
% figure;
% plot(t,SP11,'r');
% figure;
% plot(t,SP22,'g');
% figure;
% plot(t,Smx2,'r');
% figure;
% plot(t,Smx1,'g');
% figure;
% plot(t,SP11,'r');
% figure;
% plot(t,SP22,'g');

%% ------------Comparison of SDEs and EKF---------------------%%
% figure;
% plot(t,x1,'black',t,Fmx1,'red');
% figure;
% plot(t,x2,'black',t,Fmx2,'red');

% figure;
% plot(t,FP11);
% figure;
% plot(t,FP22);
%%-----------------------------------------------%%

%% ------------Comparison of EKF and Third-order Filter---------------------%%

figure;
plot(t,x1,'blue',t,Fmx1,'black',t,Tmx1,'red');
figure;
plot(t,x2,'blue',t,Fmx2,'black',t,Tmx2,'red');

figure;
subplot(2,1,1)
plot(t,TP11,'black',t,FP11,'r');
subplot(2,1,2)
plot(t,TP22,'black',t,FP22,'r');
%%-----------------------------------------------%%
