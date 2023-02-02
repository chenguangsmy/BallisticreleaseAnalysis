function [param,Fest] = sysID_2nd_Fc(dt,x,F)

d = designfilt('lowpassfir','FilterOrder',50,'PassBandFrequency',0.1,'StopBandFrequency',35,'SampleRate',1/dt);
xd = diff(x)/dt;
xd = filtfilt(d,xd);
xdd = diff(xd)/dt;

%Trim to shortest vector
x = x(3:end);
xd = xd(2:end);
F = F(3:end);

% System Id
A = [xdd' xd' x' sign(xd)'];
B = F';

param = A\B;

Fest = param(1)*xdd + param(2)*xd + param(3)*x + param(4)*sign(xd);
end