Fs=1024;
N=1024;
F0=15;
A0=0.5;
mi=2*pi;
t=[0:N-1]/Fs;
s=A0.*sin(2*pi*F0*t+mi*sin(2*pi*t));
plot(t,s);