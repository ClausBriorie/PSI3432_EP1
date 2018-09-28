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
title(s(1),'Item D: Módulo da STFD');
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
title(s(1),'Item E: Módulo da TFTD');
s(2) = subplot(212);
stem(k,angle(X_ejw));
title(s(2),'Item E: Fase da TFTD');

%% ITEM F

k=[0:1:19];

figure;

stem(k,abs(X_k));
hold on;
stem(k,abs(xk));
title('Item F: Comparação entre TFTD e SF')
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
title(s(1), 'Item G: Módulos da TDF e da TDF com um ponto a mais');
hold off;
s(2) = subplot(212);
stem(k,angle(X_k));
hold on;
stem(k2,angle(X_k2));
title(s(2), 'Item G: Fases da TDF e da TDF com um ponto a mais');
hold off;

%% ITEM H
% Mudando n para 4 períodos
n=[0:1:79];
x_n = sin(w*n*Ta) + (1/3)*sin(3*w*n*Ta);

X_k_4periodos = fft(x_n);

k=[0:1:79];

figure;
s(1) = subplot(211);
stem(k,abs(X_k_4periodos));
title(s(1),'Item H: Módulo da STFD (4 períodos)');
s(2) = subplot(212);
stem(k,angle(X_k_4periodos));
title(s(2),'Item H: Fase da SFTD (4 períodos)');

%% ITEM I
% Mudando n para 4,5 períodos
n=[0:1:89];
x_n = sin(w*n*Ta) + (1/3)*sin(3*w*n*Ta);

X_k_4emeioperiodos = fft(x_n);

k=[0:1:89];

figure;
s(1) = subplot(211);
stem(k,abs(X_k_4emeioperiodos));
title(s(1),'Item I: Módulo da STFD (4,5 períodos)');
s(2) = subplot(212);
stem(k,angle(X_k_4emeioperiodos));
title(s(2),'Item I: Fase da SFTD (4,5 períodos)');
