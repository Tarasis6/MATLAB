Fs=1024;
N=1024;
F0=80;
A0=0.5;
t=[0:N-1]/Fs;
s=A0.*sin(2*pi*F0*t);
S=fft(s);
plot(abs(S));