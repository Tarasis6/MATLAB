clear all
close all
Fs=16;
N=1024;
t=[0:N-1]/Fs;
s=sin(2*pi*t*2) + sin(2*pi*t*18+pi/2);
S=fftshift(fft(s));
ff = [-length(S)/2:length(S)/2 - 1]*Fs/N;
plot(ff,abs(S));