%exp2-PSI3432
close all;
clear all;
%EXERCICIO 1

%% ITEM A
f=500;

t=[1:1000];

x_t = sin(2*pi*f*t) + 1/3*sin(2*pi*3*f*t);

x11 = 1/(2*1i)*(exp(2*pi*f*1i));
x12 = -1/(2*1i)*(exp(-2*pi*f*1i));
x21 = 1/(6*1i)*(exp(6*pi*f*1i));
x22 = -1/(6*1i)*(exp(-6*pi*f*1i));
xk=zeros(20);

xk(2)=x11;
xk(4)=x21;
xk(18)=x22;
xk(20)=x12;

%% ITEM B


%X = 2*pi*(x11*dirac(w-2*f*pi)+x12*dirac(w+2*f*pi)+x21*dirac(w-6*f*pi) ...
%    +x22*dirac(w+6*pi*f));

%% ITEM C

Ta=0.0001;
fa=1/Ta;
f1=1/(1.9e-3);

n=[0:1:19];
w = 2*pi*f;
x_n = sin(w*n*Ta) + (1/3)*sin(3*w*n*Ta);

%% ITEM D

X_k = fft(x_n);

k=[0:1:19];

figure;
s(1) = subplot(211);
stem(k,abs(X_k));
title(s(1),'Item D: M√≥dulo da STFD');
s(2) = subplot(212);
stem(k,angle(X_k));
title(s(2),'Item D: Fase da SFTD');

%% ITEM E

X_ejw = 0;
for j=1:1:20
    X_ejw =[X_ejw x_n(j)*exp(-i*w*j)];
end

k=[0:1:20];

figure;
s(1) = subplot(211);
stem(k,abs(X_ejw));
title(s(1),'Item E: M√≥dulo da TFTD');
s(2) = subplot(212);
stem(k,angle(X_ejw));
title(s(2),'Item E: Fase da TFTD');

%% ITEM F

k=[0:1:19];

figure;
stem(k,abs(X_k));
hold on;
stem(k,abs(xk));
title('Item F: Compara√ß√£o entre TFTD e SF')
hold off;

%% ITEM G

n=[0:1:20];
w = 2*pi*f;
x_n2 = sin(w*n*Ta) + 1/3*sin(3*w*n*Ta);

X_k2 = fft(x_n2);

k=[0:1:19];
k2=[0:1:20];
figure;
s(1) = subplot(211);
stem(k,abs(X_k));
hold on;
stem(k2,abs(X_k2));
title(s(1), 'Item G: M√≥dulos da TDF e da TDF com um ponto a mais');
hold off;
s(2) = subplot(212);
stem(k,angle(X_k));
hold on;
stem(k2,angle(X_k2));
title(s(2), 'Item G: Fases da TDF e da TDF com um ponto a mais');
hold off;

%% ITEM H
% Mudando n para 4 per√≠odos
n=[0:1:79];
x_n = sin(w*n*Ta) + (1/3)*sin(3*w*n*Ta);

X_k_4periodos = fft(x_n);

k=[0:1:79];

figure;
s(1) = subplot(211);
stem(k,abs(X_k_4periodos));
title(s(1),'Item H: M√≥dulo da STFD (4 per√≠odos)');
s(2) = subplot(212);
stem(k,angle(X_k_4periodos));
title(s(2),'Item H: Fase da SFTD (4 per√≠odos)');

%% ITEM I
% Mudando n para 4,5 per√≠odos
n = [0:1:89];
x_n = sin(w*n*Ta) + (1/3)*sin(3*w*n*Ta);

X_k_4emeioperiodos = fft(x_n);

k=[0:1:89];

figure;
s(1) = subplot(211);
stem(k,abs(X_k_4emeioperiodos));
title(s(1),'Item I: M√≥dulo da STFD (4,5 per√≠odos)');
s(2) = subplot(212);
stem(k,angle(X_k_4emeioperiodos));
title(s(2),'Item I: Fase da SFTD (4,5 per√≠odos)');

%% Exerc√≠cio 2

%% ITEM C
Ta = 1/30;
fa = 1/Ta;
O_a = 2*pi*fa; % omega da amostra
N = 60;
n = [0:1:59];
res = 10000; % resolu?o

% C√°culo da TFTD te√≥rica
idx = 1;
for o = -3*pi:(pi/res):3*pi
    Y_ejomg_teo(idx) = (1/2i)*((1-exp(-1i*(o-2*pi*Ta)*N))/(1-exp(-1i*(o-2*pi*Ta))) - ...
                               (1-exp(-1i*(o+2*pi*Ta)*N))/(1-exp(-1i*(o+2*pi*Ta))));
    idx = idx + 1;
end

% C√°culo da somat√≥ria proposta
Y = zeros(1, res);
for l=1:1:3
    idx = 1;
    for o=0:(pi/res):pi - pi/res
        Y(idx) = Y(idx) + (1/2i)*(sqrt(2/pi))*((exp(1i*(o/Ta-l*O_a)+2*pi)*sinc((o/Ta-l*O_a)+2*pi)) - ...
                                               (exp(1i*(o/Ta-l*O_a)-2*pi)*sinc(2*pi-(o/Ta-l*O_a))));
        idx = idx + 1;
    end
end

Y = Y*(1/Ta);

K = 1:res;
w = -3*pi:(pi/res):3*pi;

figure;
s(1) = subplot(211);
plot(w, abs(Y_ejomg_teo));
title(s(1),'2-c): M√≥dulo da TFTD Te√≥rica');

s(2) = subplot(212);
plot(K, abs(Y));
title(s(2),'2-c): M√≥dulo da somat√≥ria proposta');

%% ex3
clear all;
close all;

L=32;

h(1:32) = 1/32;
H_k = fft(h, L);


N=1024;
H = fft(h, N);


%item a

w = 2*pi*(0:(N-1))/N;
w2 = fftshift(w);
w3 = unwrap(w2 - 2*pi);
figure;
plot(w3, abs(fftshift(H)));
xlabel('radianos');
title ('filtro H(e^j^w)');

%item b

for n=1:32
    x(n) = cos(0.1*pi*n) + 0.8*cos(0.8*pi*n);
end
figure;
plot(x);
title('x(n)');

%item c


y1 = conv(h,x);
figure;
plot (y1);

figure;
y2 = filter(H,1,x);
freqz(y2);
% plot(y2);

%Esse aqui t· dando errado sem d˙vidas

%item d

X_k = fft(x,L);
Y_k = H_k.*X_k;
y = ifft(Y_k);

figure;
plot(y);
%Esse resultado n„o faz sentido mas n„o entendi qual o problema. Comparando
%com o prÛximo resultado d· pra talvez imaginar alguma coisa (?) O que
%pensei t· no item f

%item e

L2=64;
h2(1:63) = 1/32;
H_k2 = fft(h2, L2);
for n=1:63
    x(n) = cos(0.1*pi*n) + 0.8*cos(0.8*pi*n);
end
X_k2 = fft(x,L2);
Y_k2 = H_k2.*X_k2;
y2 = ifft(Y_k2);

figure;
plot(y2);

%item f

% A diferenÁa do n˙mero de amostras È significativo para a eficiÍncia e
% funcionalidade do filtro. Com 32 coeficientes, todos os pontos do filtro
% ficaram localizados nos pontos de zero do sinc, exceto o primeiro, e 
%portanto o filtro n„o funcionou apropriadamente. Por conta disso, com mais
%coeficientes como no caso de N=63, È possÌvel verificar um sinal mais
%semelhante ao obtido utilizando a funÁ„o conv