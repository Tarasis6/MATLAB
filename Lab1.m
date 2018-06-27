A = 1.2;
tau = 5;
T = 3*tau;
t = 0:0.01:T;
f = -50:0.01:50;
w = 2*pi*f; 
s = A.*(t<tau); % прямоугольный импульс
A1 = abs(s);
% спектр прямоугольного сигнала (определён как функция)
Upuls=inline('tau*sin(w*tau/2)./(w*tau/2).*exp(-32i*w*tau/2)',...
 'w','tau');
% вычисляем спектры сигналов с помощью свойств
S = Upuls(w, tau); 
%вычисление спектра дискретного сигнала
X2 = zeros(1, length(f));
for n=1:length(t)
 X2 = X2 + 0.1*s(n)*exp(-i*w*(n-1)*0.02);
end
S1 = fftshift(fft(s));
%сравниваем амплитудные спектры
subplot(3,2,1), plot(f,abs(S)), xlabel('f'), ylabel('anal spectr')
subplot(3,2,2), stem(f,abs(X2)), xlabel('f'), ylabel('diskr spectr')
subplot(3,2,3), plot(t,s), xlabel('t'),ylabel('A')
subplot(3,2,6), plot (f,abs(S)), xlabel('f')
subplot(3,2,5), plot(abs(S1)), xlabel('fft')
subplot(3,2,4), stem(t,A1)
