Fs=16;
N=1024;
t1=[0:31*16-1]/16/Fs;
t=[0:31]/Fs;
s=sin(2*pi*t*2);
s1=sin(2*pi*t1*(32+2));
plot(t1,s1,'r');
hold on
plot(t,s)
plot(t,s,'.')
hold off