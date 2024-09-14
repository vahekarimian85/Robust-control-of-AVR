%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%                      Robust control of AVR                          %%%
%%%                        (H2/Hinf method)                             %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all

s=tf('s')
%% define weighting matrices
a=1;b=0.3;c=500;d=0.3;
e=1;f=4.5;g=3;h=2.8;  
Ws=(a*s+b)/(c*s+d);
Wu=(e*s+f)/(g*s+h);
Ws_1=1/Ws;
figure
bodemag(Ws_1)
legend('1/Wp')
title('Bode mag. diagram of 1/Wp')
figure
bodemag(Wu)
legend('Wu')
title('Bode mag.diagram of Wu')

%% interconnected system for nominal plant
Wn=1/(0.0305*s+1);
G0=(116.8/(0.06*s+1))*(1/(0.7*s+1))*(1/(1.5*s+1));
systemnames='G0 Ws Wu Wn';
inputvar='[r;d;n;u]';
outputvar='[Ws;Wu;r-Wn]';
input_to_G0='[u]';
input_to_Ws='[r-Wn]';
input_to_Wn='[n+G0]';
input_to_Wu='[u]';
sysoutname='P_ic';
cleanupsysic='yes';
sysic
P_ic=minreal(P_ic);

%% design H-inf controller
wz=[1 1];
[k,cl,gg,info]=h2hinfsyn(ss(P_ic),1,1,1,wz);
tf(k)
zpk(k)
figure
sigma(cl)
title('Singular values of closed loop system')
figure
bodemag(k)
title('Bode mag. of H2/Hinf controller')

%% reference input step response
Wn_1=1/(0.006*s+1);
G1=(116.8/(0.06*s+1))*(1/(0.7*s+1))*(1/(1.5*s+1));
GG=k*G1/(1+G1*k*Wn_1);
figure
step(G1)
title('Step response of open loop sys.')
figure
step(GG)
title('Step response of closed loop sys.')
figure
subplot(2,1,1)
bodemag(G1)
title('Bode mag. of open loop sys.')
subplot(2,1,2)
bodemag(GG)
title('Bode mag. of closed loop sys.')

%% stste space and transfer function representation of interconnected sys.
TA=0.06;TE=0.7;TG=1.5;KS=1;TS=0.006;PTA=0.67;PTE=0.43;PTG=1/3;KF=203.5;PKF=0.966;
z1=a/c;z2=(b/d)-(a/c);z3=c/d;x1=e/g;x2=(f/h)-(e/g);x3=g/h;
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

%% Mu plot (obtained final Mu)
block=[1 0;1 0;1 0;1 0;3 2];
[bounds,muinfo]=mussv(Tzw,block);
bounds
figure
semilogx(bounds(:,2))
title('Mu plot of closed loop system with H2/Hinf controller')

%% time responses
figure
subplot(2,1,1)
simout=sim('AVR_H2_Hinf__')
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
simout=sim('AVR_H2_Hinf_nonlinear_')
plot(simout.vv(:,1))
xlabel('time(sec)');ylabel('V(t)')
title('Output Voltage')
hold on
subplot(2,1,2)
plot(simout.uu(:,1))
xlabel('time(sec)');ylabel('U(t)')
title('control input')
