%first order
j = sqrt(-1);
w = 0:2*pi*0.01:2*pi*2;

a0 = 1;
a1 = -0.8;
b0 = 1;
b1 = -1.25;
A = [b0 a1];
B = [a0 b1];

%H1 = b0+b1*z^(-1)/1+a1*z^(-1);

%Frequency response

figure (1);

H1=(b0-b1.*exp(-j.*w))./(1+a1.*exp(-j.*w));
h1=abs(H1);
subplot(2,2,1);
plot(w,h1);
xlabel('w, рад');ylabel('АЧХ');

%Phase response

f1 = angle(H1);
subplot(2,2,2)
plot(w,f1);
xlabel('w, рад');ylabel('ФЧХ');

%Impulse response

[h1,n1]=impz(B,A); 
subplot(2,2,3);
stem(n1,h1);
xlabel('n');ylabel('h(n)');

%zeros and poles

num = [b0 b1];
den = [a0 a1];
subplot(2,2,4);
zplane(num, den); 

%second order

a0 = 1;
a1 = 0;
a2 = 0;
b0 = 0.5;
b1 = 1;
b2 = 0.5;

%H1 = b0+b1*z^(-1)+b2*z^(-2)/1+a1*z^(-1)+a2*z^(-2);

%Frequency response

figure (1);

H1=(b0+b1.*exp(-j.*w))./(1+a1.*exp(-j.*w));
h1=abs(H1);
subplot(2,2,1);
plot(w,h1);
xlabel('w, рад');ylabel('АЧХ');

%Phase response

f1 = angle(H1);
subplot(2,2,2)
plot(w,f1);
xlabel('w, рад');ylabel('ФЧХ');

%Impulse response

[h1,n1]=impz(B,A); 
subplot(2,2,3);
stem(n1,h1);
xlabel('n');ylabel('h(n)');


num = [b0 b1 b2];
den = [a0 a1 a2];
zplane(num, den) 

