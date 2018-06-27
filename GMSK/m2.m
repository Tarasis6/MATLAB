beta = 0.3;
Fs = 1000;
Br = 200;
Tr = 1/Br;
K = 1;
F0 = 8*Br;
T = 16 ;

t = -K*Tr:1/Fs:K*Tr;
g = beta/Tr*sqrt(2*pi/log(2))*exp(-2*(pi*beta)^2*t.^2/log(2)/(Tr^2));
g = g/sum(g);
data = randi(2,1,T*Br)-1;
data = reshape(repmat(data,Fs/Br,1),[],1)*2 - 1;

data_f = conv(data,g);
data_f = data_f((length(g)-1)/2+1:end-(length(g)-1)/2);
plot(data_f)
hold on
plot(data,'r')
hold off
figure

t = [0:T*Fs-1]'/Fs;


wd= pi*Br*2;

smi=cumsum(data);
integral=smi*1/Fs;
integral=integral*wd;
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

data  = data_f;
t = [0:T*Fs-1]'/Fs;


wd= pi*Br*2;

smi=cumsum(data);
integral=smi*1/Fs;
integral=integral*wd;
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