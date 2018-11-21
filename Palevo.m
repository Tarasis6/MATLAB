close all;
clear all

N = 2^13;
m  = 2;
Br = 16;
Fs = 256;
f0 = 32;
beta = 0.7;
t = [0 : Fs / Br * N - 1]' / Fs;
%%%%%%%%%%%%5
dF = 0.3;
fi = 0 ;

dFk = 0.0;
fik = 0;


data = randi(2,N,m)-1;
table_qpsk = [1 + 1i,1 - 1i,-1 + 1i, -1 - 1i];
IQ  = table_qpsk(bi2de(data) + 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IQf = [IQ;zeros(Fs/Br-1,length(IQ))];
IQf = reshape(IQf,[],1);
IQ = repmat(IQ,Fs/Br,1);
IQ = reshape(IQ,[],1);

b = rcosdesign(beta,4,Fs/Br,'sqrt') / 4;
IQf = conv(IQf,b,'same');
Sf = real(IQf) .* cos(2 * pi * f0 * t) - imag(IQf) .* sin(2 * pi * f0 * t);

S = real(IQ) .* cos(2 * pi * f0 * t) - imag(IQ) .* sin(2 * pi * f0 * t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure
% subplot(2,1,1)
% plot(t(1:Fs/Br*8-1),S(1:Fs/ Br*8-1));
% hold on;
% plot(t(1:Fs/Br*8-1),real(IQ(1:Fs/Br*8-1)),'r');
% plot(t(1:Fs/Br*8-1),imag(IQ(1:Fs/Br*8-1)),'g');
% hold off;
% 
% subplot(2,1,2)
% plot(t(1:Fs/Br*8-1),Sf(1:Fs/Br*8-1));
% hold on;
% plot(t(1:Fs/Br*8-1),real(IQf(1:Fs/Br*8-1)),'r');
% plot(t(1:Fs/Br*8-1),imag(IQf(1:Fs/Br*8-1)),'g');
% hold off;
% 
% figure
% subplot(1,2,1)
% plot(real(IQ),imag(IQ));
% axis([-2,2,-2,2]);
% subplot(1,2,2)
% plot(real(IQf*Fs/Br),imag(IQf*Fs/Br));
% axis([-2,2,-2,2]);
% 
% figure
% %***********************************************************************
% ff = [-length(S) / 2 : length(S) / 2 - 1 ] * Fs / length(S); %підпис осі
% %***********************************************************************
% Spec = fftshift(fft(IQ));
% Specf = fftshift(fft(IQf));
% subplot(1,2,1)
% plot(ff,abs(Spec));
% subplot(1,2,2)
% plot(ff,abs(Specf));
% 
% figure
% Spec = fftshift(fft(S));
% Specf = fftshift(fft(Sf));
% subplot(1,2,1)
% plot(ff,abs(Spec));
% subplot(1,2,2)
% plot(ff,abs(Specf));
% %=======================================================================
% 
% 
% figure
% plot(real(IQf*Fs/Br),'r');
% hold on;
% plot(imag(IQf*Fs/Br),'g');
% 
% plot(real(IQ),'r');
% plot(imag(IQ),'g');
% 
% axis([0,150,-1.5,1.5]);
% 
% hold off;
% 
% %*************** GLAS diagram   *************************
% IQ_eye = reshape(IQf,4 * Fs/Br,[]);
% figure
% plot(real(IQ_eye));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IQr = IQ.*exp(1i*(2*pi*dF*t - fi));
IQrf = IQf.*exp(1i*(2*pi*dF*t - fi));

figure
%subplot(1,2,1)
plot(real(IQrf(1:Fs/Br:end)*Fs/Br),imag(IQrf(1:Fs/Br:end)*Fs/Br),'.');
%axis([-2,2,-2,2]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% IQr = IQr.*exp(-1i*(2*pi*dFk*t - fik));
%IQrf = IQrf.*exp(-1i*(2*pi*dFk*t - fik));


% subplot(1,2,2)
% plot(real(IQrf(1:Fs/Br:end)*Fs/Br),imag(IQrf(1:Fs/Br:end)*Fs/Br),'.');
% axis([-2,2,-2,2]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IQrr = (IQrf(1:Fs/Br:end)).^4;


b = rcosdesign(beta,4,Fs/Br,'sqrt') / 4;
IQf = conv(IQrf,b,'same');

hFIR = dsp.FIRFilter
hFIR.Numerator  = rcosdesign(beta,4,Fs/Br,'sqrt') / 4;
n = 0;
k=1;
while(n < length(IQrf)/Fs/Br)
    
    IQrf2 = step(hFIR, IQrf(n*FS/Br+1 : (n+1)*Fs*Br));
    
    if((imag(IQrf2(k))>0) && (real(IQrf2(k))>0))
    alfa(n) = atan(imag(IQrf2(k))/real(IQrf2(k))) - pi/4 + 2*pi;
    end
    if((imag(IQrf2(k))<0) && (real(IQrf2(k))<0))
    alfa(n) = atan(imag(IQrf2(Fs/Br))/real(IQrf2(Fs/Br))) - 3*pi/4 + 2*pi;
    end
    if((imag(IQrf(k))2>0) && (real(IQrf2(k))<0))
    alfa(n) = atan(imag(IQrf2(Fs/Br))/real(IQrf2(Fs/Br))) - 5*pi/4 + 2*pi;
    end
    if((imag(IQrf2(k))<0) && (real(IQrf2(k))>0))
    alfa(n) = atan(imag(IQrf2(Fs/Br))/real(IQrf2(Fs/Br))) - 3*pi/4 + 2*pi;
    end
end
figure
plot(real(IQrr),imag(IQrr),'.');
%alfa = (IQrr(1:end-1)'*IQrr(2:end))/(length(IQrr)-1);
%angle(alfa)/4*Br/2/pi
alfa = angle(atan(imag(IQrr)./real(IQrr)));
figure
%plot(real(alfa),imag(alfa),'.');
%axis([-6,6,-6,6]);


%***********************************************************************

SSS = fftshift(fft(IQrr));
ff = [-length(SSS) / 2 : length(SSS) / 2 - 1 ] * Fs / length(SSS); %підпис осі
figure
plot(ff,abs(SSS));
[v,c] =  max(abs(SSS));
dFc = ff(c)/4


IQrf = IQf.*exp(1i*(2*pi*dF*t - fi));
IQrf = IQrf.*exp(-1i*(2*pi*dFc*t - fik));

figure
plot(real(IQrf(1:Fs/Br:end)*Fs/Br),imag(IQrf(1:Fs/Br:end)*Fs/Br),'.');
%axis([-2,2,-2,2]);




