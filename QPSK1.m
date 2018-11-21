close all;
clear all

N = 2^12;
m  = 2;
Br = 16;
Fs = 32;
f0 = 32;
beta = 0.5;
t = [0 : Fs / Br * N - 1]' / Fs;
dF = 0.1;
fi = 0

PhaseRecoveryLoopBandwidth = 0.05;
PostFilterOversampling = 2; 

data = randi(2,N,m)-1;
table_qpsk = [1 + 1i,1 - 1i,-1 + 1i, -1 - 1i]/sqrt(2);
IQ  = table_qpsk(bi2de(data) + 1);

hSRCFilt = dsp.FIRFilter('Numerator',rcosdesign(beta,8,Fs/Br,'sqrt'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IQf = [IQ;zeros(Fs/Br-1,length(IQ))];
IQf = reshape(IQf,[],1);
IQ = repmat(IQ,Fs/Br,1);
IQ = reshape(IQ,[],1);

b = rcosdesign(beta,8,Fs/Br,'sqrt');
IQf = conv(IQf,b,'same')/4;
% % % % Sf = real(IQf) .* cos(2 * pi * f0 * t) - imag(IQf) .* sin(2 * pi * f0 * t);
% % % % 
% % % % S = real(IQ) .* cos(2 * pi * f0 * t) - imag(IQ) .* sin(2 * pi * f0 * t);
% % % % 
% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % figure
% % % % subplot(2,1,1)
% % % % plot(t(1:Fs/Br*8-1),S(1:Fs/Br*8-1));
% % % % hold on;
% % % % plot(t(1:Fs/Br*8-1),real(IQ(1:Fs/Br*8-1)),'r');
% % % % plot(t(1:Fs/Br*8-1),imag(IQ(1:Fs/Br*8-1)),'g');
% % % % hold off;
% % % % 
% % % % subplot(2,1,2)
% % % % plot(t(1:Fs/Br*8-1),Sf(1:Fs/Br*8-1));
% % % % hold on;
% % % % plot(t(1:Fs/Br*8-1),real(IQf(1:Fs/Br*8-1)),'r');
% % % % plot(t(1:Fs/Br*8-1),imag(IQf(1:Fs/Br*8-1)),'g');
% % % % hold off;

figure
subplot(1,2,1)
plot(real(IQ),imag(IQ));
axis([-2,2,-2,2]);
subplot(1,2,2)
plot(real(IQf*Fs/Br),imag(IQf*Fs/Br));
axis([-2,2,-2,2]);

figure
%***********************************************************************
ff = [-length(IQ) / 2 : length(IQ) / 2 - 1 ] * Fs / length(IQ); %підпис осі
%***********************************************************************
Spec = fftshift(fft(IQ));
Specf = fftshift(fft(IQf));
subplot(1,2,1)
plot(ff,abs(Spec));
subplot(1,2,2)
plot(ff,abs(Specf));

% % % figure
% % % Spec = fftshift(fft(S));
% % % Specf = fftshift(fft(Sf));
% % % subplot(1,2,1)
% % % plot(ff,abs(Spec));
% % % subplot(1,2,2)
% % % plot(ff,abs(Specf));
%=======================================================================


figure
plot(real(IQf*Fs/Br),'r');
hold on;
plot(imag(IQf*Fs/Br),'g');

plot(real(IQ),'r');
plot(imag(IQ),'g');

axis([0,150,-1.5,1.5]);

hold off;

%*************** GLAS diagram   *************************
IQ_eye = reshape(IQf,4 * Fs/Br,[]);
figure
plot(real(IQ_eye));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IQr = IQ.*exp(1i*(2*pi*dF*t - fi));
IQrf = IQf.*exp(1i*(2*pi*dF*t - fi));

figure
subplot(1,2,1)
plot(real(IQr(1:Fs/Br:end)),imag(IQr(1:Fs/Br:end)),'.');
axis([-2,2,-2,2]);
subplot(1,2,2)
plot(real(IQrf(1:Fs/Br:end)*Fs/Br),imag(IQrf(1:Fs/Br:end)*Fs/Br),'.');
axis([-2,2,-2,2]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PhaseDDCOut = 0;
j=1;
Integration = 0;
loopFilter = 0;

PhaseRecoveryDampingFactor = 1/sqrt(2);



K = 1;
A = 1/sqrt(2);
PhaseErrorDetectorGain = 2*K*A^2+2*K*A^2; % K_p for Fine Frequency Compensation PLL, determined by 2KA^2 (for binary PAM), QPSK could be treated as two individual binary PAM
PhaseRecoveryGain = 1; % K_0 for Fine Frequency Compensation PLL

theta = PhaseRecoveryLoopBandwidth/...
    (PhaseRecoveryDampingFactor + ...
    0.25/PhaseRecoveryDampingFactor)/PostFilterOversampling;
d = 1 + 2*PhaseRecoveryDampingFactor*theta + theta*theta;
Kproportional = (4*PhaseRecoveryDampingFactor*theta/d)/...
    (PhaseErrorDetectorGain*PhaseRecoveryGain);
Kintegral = (4*theta*theta/d)/...
    (PhaseErrorDetectorGain*PhaseRecoveryGain);

figure
SR = zeros(100,1);
for i=1:length(IQrf)
    IQrf_fc(i) = IQrf(i)*exp(-1i*2*pi*PhaseDDCOut(i));
    IQrf_filt(i) = step(hSRCFilt,IQrf_fc(i));
    
    if (mod(i-1,Fs/Br) == 0)
        e = imag(IQrf_filt(i))*sign(real(IQrf_filt(i))) - real(IQrf_filt(i))*sign(imag(IQrf_filt(i)));
        
        DATA_OUT = IQrf_filt(i);
        j=j+1;
    else
        e=0;
    end
    PhaseDDCOut(i+1) = PhaseDDCOut(i) + loopFilter;
    Integration = Integration + e*Kintegral;
    loopFilter = Integration + e*Kproportional;
    
    SR = [SR(2:end); DATA_OUT];
    plot(real(SR),imag(SR),'.');
    axis([-0.5,0.5,-0.5,0.5])
    drawnow
    
end



