%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%                   Parameteric controller for AVR                    %%%
%%%                                                                     %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all

s=tf('s')

%% parameteric controller
k=(0.063*s^4+1.245*s^3+3.442*s^2+3.26*s+1)/(0.063*s^4+8.54*s^3+140.3*s^2+148.2*s);
k1=tf(k)
omega=logspace(-3,4,200);
k2=frd(k1,omega);

%% bode mag. of controller and system
G=116.8/((0.06*s+1)*(0.7*s+1)*(1.5*s+1));
Wn=1/(0.006*s+1);
Wn_1=1/Wn;
f=G*k/(1+G*k*Wn_1);
figure
bodemag(k)
title('Bode mag. of parameteric controller')
figure
subplot(2,1,1)
bodemag(G)
title('Bode mag. of open loop system')
subplot(2,1,2)
bodemag(f)
title('Bode mag. of closed loop system')

%% stste space and transfer function representation of interconnected sys.
a=1;b=1;c=1;d=1;e=1;f=1;g=1;h=1;
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
Tzw=M11+M12*k2*inv(eye(1)-M22*k2)*M21;

%% Mu plot (obtained final Mu)
block=[1 0;1 0;1 0;1 0;3 2];
[bounds,muinfo]=mussv(Tzw,block);
bounds
figure
semilogx(bounds(:,2))
title('Mu plot of closed loop system with parameteric controller')

%% time responses
figure
subplot(2,1,1)
simout=sim('AVR_parameteric__')
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
simout=sim('AVR_parameteric_nonlinear_')
plot(simout.vv(:,1))
xlabel('time(sec)');ylabel('V(t)')
title('Output Voltage')
hold on
subplot(2,1,2)
plot(simout.uu(:,1))
xlabel('time(sec)');ylabel('U(t)')
title('control input')
