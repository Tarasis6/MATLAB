Fs = 1024;
Br = 16;
F0 = 8*Br;
T = 8;
data = randi(1,T*Br);
data = reshape(repmat(data,Fs/Br,1),[],1)*2 - 1;
t = [0:T*Fs-1]'/Fs;

wd= pi*Br*1/2;

smi=cumsum(data);
integral=smi*1/Fs;
integral=integral.*wd;
I=cos(integral);
Q=sin(integral);
s=I.*cos(2*pi*F0*t)-Q.*sin(2*pi*F0*t);
figure
plot(t,s);
hold on
plot(t,data,'r')
hold off


figure
ff = [-length(s)/2:length(s)/2-1]/T;
semilogy(ff,abs(fftshift(fft(s))))
