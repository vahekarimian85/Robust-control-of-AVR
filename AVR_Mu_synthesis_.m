%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%                      ROBUST CONTROL OF AVR                          %%%
%%%                         (Mu Synthesis)                              %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all

s=tf('s')
%% define weighting matrices
a=1;b=220;c=210000;dd=1;  
ee=1;f=1;g=1;h=20;       
Ws=(a*s+b)/(c*s+dd);    
Wu=(ee*s+f)/(g*s+h);
Ws_1=1/Ws;
omega=logspace(-2,3,200);
figure
bodemag(Ws_1)
legend('1/Wp')
figure
bodemag(Wu,omega)
legend('Wu')
Wn=1/(0.006*s+1);

%% constructing interconnected system for uncertain plant
G0=(203.5/(0.06*s+1))*(1/(0.7*s+1))*(1/(1.5*s+1));
delta=ultidyn('delta',[1 1]);
W=10*(s+1)/(s+10);
G=G0*(1+W*delta);
M=iconnect;
r=icsignal(1);
d=icsignal(1);
n=icsignal(1);
e=icsignal(1);
u=icsignal(1);
zs=icsignal(1);
zu=icsignal(1);
M.Equation{1}=equate(e,r-Wn*(n+d+G*u));
M.Input=[r;d;n;u];
M.Output=[Ws*e;Wu*u;e];
T=tf(M.System)

%% designing controller by mu-synthesis method
[k,cl,bnd,dkinfo]=dksyn(M.System,1,1)
tf(k)
zpk(k)
figure
bodemag(k)
legend('K')
title('Bode mag. of Mu-synthesis controller')

%% reference input step response
G1=(203.5/(0.06*s+1))*(1/(0.7*s+1))*(1/(1.5*s+1));
GG=k*G1/(1+G1*k*Wn);
figure
step(G1)
title('step response of open loop sys.')
figure
step(GG)
title('step response of closed loop sys.')
figure
subplot(2,1,1)
bodemag(G1)
title('Bode mag. of open loop sys.')
subplot(2,1,2)
bodemag(GG)
title('Bode mag. of closed loop sys.')

%% mu-plot of each iteration
figure
l1=semilogx(dkinfo{1,1}.MussvBnds(:,2),'k--')
hold on
l1.LineWidth=2;
l2=semilogx(dkinfo{1,2}.MussvBnds(:,2),'g-.')
l2.LineWidth=2;
hold on
l3=semilogx(dkinfo{1,3}.MussvBnds(:,2),'r:')
l3.LineWidth=2;
legend('\mu_1','\mu_2','\mu_3')

%% state-space and transfer function representation of interconnected sys.
TA=0.06;TE=0.7;TG=1.5;KS=1;TS=0.006;PTA=0.67;PTE=0.43;PTG=1/3;KF=203.5;PKF=0.966;
z1=a/c;z2=(b/dd)-(a/c);z3=c/dd;x1=ee/g;x2=(f/h)-(ee/g);x3=g/h;
A=[-1/TA 0 0 0 0 0
    1/TE -1/TE 0 0 0 0
    0 1/TG -1/TG 0 0 0
    0 0 KS/TS -1/TS 0 0
    0 0 0 -z2/z3 -1/z3 0
    0 0 0 0 0 -1/x3];
B1=[1/TA -PTA/TA 0 0 0 0 0
    0 0 -PTE/TE 0 0 0 0
    0 0 0 -PTG/TE 0 0 0
    0 0 0 0 0 KS/TS KS/TS 
    0 0 0 0 z2/z3 0 0
    0 0 0 0 0 0 0];
C1=[0 0 0 0 0 0
    -1 0 0 0 0 0
    1 -1 0 0 0 0
    0 1 -1 0 0 0
    0 0 0 -z1 1 0
    0 0 0 0 0 1];
B2=[KF/TA;0;0;0;0;x2/x3];
C2=[0 0 0 -1 0 0];
D11=[0 0 0 0 0 0 0
    1 -PTA 0 0 0 0 0
    0 0 -PTE 0 0 0 0
    0 0 0 -PTG 0 0 0
    0 0 0 0 z1 0 0
    0 0 0 0 0 0 0];
D12=[KF*PKF;KF;0;0;0;x1];
D21=[0 0 0 0 1 0 0];
D22=[0];
B=[B1 B2];
C=[C1;C2];
D=[D11 D12;D21 D22];
M2=C*inv(s*eye(6)-A)*B+D;
M11=M2(1:6,1:7);
M12=M2(1:6,8);
M21=M2(7,1:7);
M22=M2(7,8);

%% lower LFT ( F_l(M,K) )
Tzw=M11+M12*k*inv(eye(1)-M22*k)*M21;

%% mu-plot (obtained final mu)
block=[1 0;1 0;1 0;1 0;3 2];
[bounds,muinfo]=mussv(Tzw,block);
bounds
figure
semilogx(bounds(:,2))
axis([0.00001 10000 0 3])

%% time responses
figure
subplot(2,1,1)
simout=sim('AVR_Mu_synthesis__')
plot(simout.v(:,1))
xlabel('time(sec)');ylabel('V(t)')
title('Output Voltage')
hold on
subplot(2,1,2)
plot(simout.u(:,1))
xlabel('time(sec)');ylabel('U(t)')
title('control input')

%% time response of real nonlinear system
figure
subplot(2,1,1)
simout=sim('AVR_Mu_nonlinear_')
plot(simout.vv(:,1))
xlabel('time(sec)');ylabel('V(t)')
title('Terminal Voltage')
hold on
subplot(2,1,2)
plot(simout.uu(:,1))
xlabel('time(sec)');ylabel('U(t)')
title('control input')