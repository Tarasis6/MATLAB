Fs=1024;
N=1024;
F0=80;
mi=0.5;
t=[0:N-1]/Fs;
s=(1+mi*sin(2*pi*8*t)).*sin(2*pi*F0*t);
S=fft(s);
plot(fftshift(abs(S)));