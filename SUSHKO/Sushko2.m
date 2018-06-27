close all; 
clear all;
 

Z=-1000:2000;
w=0:2000;
 
nun1=[-0.035 -0.105 -0.105 -0.035 ];
den1=[-1 1.4 -0.863 0.184536];
 
[h2,n2]=impz(nun1,den1);
 
K = 6.483e8./(((i*w)+865.477).*((i*w).^2+865.477*(i*w)+7.491e5));
plot(w,abs(K));
grid on
 
Z=exp((-i.*w*2*pi)/1000);
H2=(-0.035-0.105*Z.^(1)-0.105*Z.^(2)-0.035*Z.^(3))./(-1+1.4*Z.^(1)-0.863*Z.^(2)+0.184536*Z.^(3));
 
figure(1);
subplot(2,2,1);
zplane(nun1,den1);
grid on
 
subplot(2,2,2);
stem(n2,h2);
grid on
 
figure(1);
subplot(2,2,4);
plot(w,abs(K));
grid on
 
subplot(2,2,3);
plot(w,abs(H2));
grid on
 
 
 
 

