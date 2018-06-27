Fs=1e6;
N=1024;
t=[0:N-1]/Fs;
S=3*sin(2*pi*1000*t);
plot(t,S);