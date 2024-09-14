clc;
clear
close all

s=tf('s')
G0=(210.5/(0.06*s+1))*(1/(0.7*s+1))*(1/(1.5*s+1));
figure
for KF=7:196.5:400
    for TA=0.02:0.04:0.1
        for TE=0.4:0.3:1
            for TG=1:1:2
                G=(KF/(TA*s+1))*(1/(TE*s+1))*(1/(TG*s+1));
                W=(G/G0)-1;
                bodemag(W)
                hold on
            end
        end
    end
end
hold on
bodemag(15*(s+1)/(s+10),'*r')

