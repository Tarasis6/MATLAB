Fs=16;
N=1024;
t=[0:N-1]/Fs;
s=sin(2*pi*t*2);
S=fftshift(fft(s));
plot(abs(S));