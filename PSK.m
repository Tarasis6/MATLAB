close all
T = 10  ;
Sr = 16;
Fs = 1024;
F0 = 128;
n = 3; 
t = [0:Fs*T-1]'/Fs;
N = Fs*T;
Br = Sr*n;
data = randi(2,T*Sr,n)-1;
IQ_table = [1+0*i, 1/sqrt(2)+1/sqrt(2)*i, 0+1*i, -1/sqrt(2)+1/sqrt(2)*i, -1+0*i, -1/sqrt(2)-1/sqrt(2)*i, 0-1*i, 1/sqrt(2)-1/sqrt(2)*i ];
data = bi2de(data) + 1;
symbol = IQ_table(data);

data = reshape(repmat(symbol,Fs/Sr,1),[],1);
s = real(data).*cos(2*pi*F0*t)-imag(data).*sin(2*pi*F0*t);
I = real(data);
Q = imag(data);
plot(t,s);
hold on
plot(t,I,'r');
plot(t,Q,'g');
hold off
figure
plot(I,Q);
figure
eyeI = reshape(I,T,[]);
eyeQ = reshape(Q,T,[]);
plot(eyeI,'b')
figure
ff = [-N/2:N/2-1]*Fs/N;
plot(ff,abs(fftshift(fft(s))))

hFilt = rcosine(Sr,Fs,'sqrt',0.5);
Nfilt = length(hFilt);
    
Ifilt = conv(I,hFilt);
Qfilt = conv(Q,hFilt);
Ifilt = Ifilt((Nfilt-1)/2+1:end - (Nfilt-1)/2);
Qfilt = Qfilt((Nfilt-1)/2+1:end - (Nfilt-1)/2);

sFilt = Ifilt.*cos(2*pi*F0*t)-Qfilt.*sin(2*pi*F0*t);
figure
plot(t,s)
hold on
plot(t,Ifilt,'r');
plot(t,Qfilt,'g')
hold off

figure
plot(t,I*80,'r')
hold on
plot(t,Ifilt,'r');
plot(t,Q*80,'g')
plot(t,Qfilt,'g')
hold off
figure
plot(Ifilt,Qfilt,'.');

figure
plot(Ifilt(32:Fs/Sr:end),Qfilt(32:Fs/Sr:end),'.');

figure
ff = [-N/2:N/2-1]*Fs/N;
plot(ff,abs(fftshift(fft(sFilt/10))))
hold on
plot(ff,abs(fftshift(fft(s))),'r')
hold off
