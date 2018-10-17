%% Exp - 3 - PSI3432
close all;
clear all;
% Exerc?cio 1 - Decima??o
%% Item 1 - Filtro usando janelas de Kaiser

delta_p = 0.02;
delta_r = 0.01;
delta_omega = (13e3 - 11e3)*((2*pi)/24e3);

p = -20*log10(delta_p);
r = -20*log10(delta_r);
A = max(p, r);

beta = 0.5842*(A-21)^(0.4) + 0.07886*(A-21); % Pois A=40 (21<A<50)

N = floor(((A-8)/(2.285*delta_omega))+1);
janela_kaiser = kaiser(N, beta);
figure;
plot(janela_kaiser)
title('Ex1-1: Janela de Kaiser projetada');

%% Item 2 - Sinal x[n] ap?s passar pelo filtro
fa = 48; % [kHz]
Ta = 1/fa;
duracao = 3; % [ms]

t_am = 0:Ta:duracao;
x_n= sin(2*pi*t_am) - (1/3)*sin(42*pi*t_am);

figure;
s(1) = subplot(211);
plot(t_am, x_n);
title(s(1),'Ex1-2: Sinal x[n]');

G = 12.30; % Ajuste de ganho do filtro
y_n = filter(janela_kaiser, G, x_n);

s(2) = subplot(212);
plot(t_am, y_n);
title(s(2),'Ex1-2: Sinal apos filtragem com janela Kaiser');


%% Item 3 - Compara??es dos ganhos

% Reconstru??o do sinal, mas com dura??o aumentada
fa = 48;
Ta = 1/fa;
duracao = 13;

t_am = 0:Ta:duracao;
x_n = sin(2*pi*t_am) - (1/3)*sin(42*pi*t_am);

G = 12.30; % Ajuste de ganho do filtro
y_n = filter(janela_kaiser, G, x_n);


% Elimina??o do transit?rio do filtro
L = 38; % Atraso do filtro
y_n = y_n(2*L:2*L+480-1);
t_am_y = t_am(2*L:2*L+480-1);

% C?lculos das transformadas
x_n = x_n(1:480);
X_ejw = fft(x_n);

Y_ejw = fft(y_n);

figure;
s(1) = subplot(211);
stem(abs(X_ejw));
title(s(1),'Ex1-3: Modulo da TDF de x');

s(2) = subplot(212);
stem(abs(Y_ejw));
title(s(2),'Ex1-3: Modulo da TDF de y');

% A partir da an?lise dos gr?ficos, percebe-se que a parcela do sinal que
% apresenta frequ?ncia mais alta, de 21 kHz, ? praticamente eliminada no
% no sinal de sa?da

%% Item 4 - Redu??o da taxa de amostragem

fa = 48;
Ta = 1/fa;
duracao = 13;

t_am = 0:Ta:duracao;
x_n = sin(2*pi*t_am) - (1/3)*sin(42*pi*t_am);

G = 12.30;
y_n = filter(janela_kaiser, G, x_n);

y_n_24kHz = zeros(1, 312);
for i=1:312
    y_n_24kHz(i) = y_n(2*i);
end

figure;
s(1) = subplot(211);
plot(t_am, x_n);
title(s(1),'Ex1-4: Sinal com taxa de amostragem reduzida');

s(2) = subplot(212);
plot(y_n_24kHz);
title(s(2),'Ex1-4: Saida do filtro x[n]');


% Exerc?cio 2 - Interpola??o
%% Item 1 - Filtro

delta_p = 0.015;
delta_r = 0.01;
delta_omega = 2e3*2*pi/24e3;

p = -20*log10(delta_p);
r = -20*log10(delta_r);
A = max(p, 40);

beta = 0.5842*(A-21)^(0.4) + 0.07886*(A-21); % Pois A=40 (21<A<50)

N = floor(((A-8)/(2.285*delta_omega))+1);
janela_kaiser = kaiser(N, beta);
figure;
plot(janela_kaiser)
title('Ex2-1: Janela de Kaiser projetada');


%% Item 2
fa = 48; % [kHz]
Ta = 1/fa;
duracao = 3; % [ms]

t_am = 0:Ta:duracao;
x_n= sin(2*pi*t_am) - (1/3)*sin(42*pi*t_am);
y_l = zeros(1, 312);
for n=1:1:141
    y_l(3*n)=x_n(n);
    y_l(3*n+1)=0;
    y_l(3*n+2)=0;
end

G=1.8; %calculado experimentalmente para gerar onda com ganho 3

xlinha_n = filter(janela_kaiser, G, y_l);
figure;
plot(xlinha_n);
title('Ex2-2: xlinha');

%% Item 3
fa = 144; % [kHz]
Ta = 1/fa;
duracao = 3; % [ms]

t_am = 0:Ta:duracao;
x_n2= sin(2*pi*t_am) - (1/3)*sin(42*pi*t_am);

figure;
plot(x_n2)
hold on;
plot(xlinha_n);
hold off;
legend('ideal', 'experimental');
title('Ex2-3: Comparacao de xlinha e o sinal ideal');


%% Item 4
X_ejw2 = fft(xlinha_n);
Y_ejw2 = fft(janela_kaiser);

figure;
s(1) = subplot(211);
stem(abs(X_ejw2));
title(s(1),'Ex1-3: Modulo da TDF de xlinha');

s(2) = subplot(212);
stem(abs(Y_ejw2));
title(s(2),'Ex1-3: Modulo da TDF do filtro');

%% Resposta da observação 1:
%    As operações de interpolação e decimação não são comutativas
%    em geral, pois a decimação apresenta potencial de destruição de
%    informação contida no sinal. Logo, exceto em casos nos quais os
%    fatores de interpolação e decimação sejam co-primos, é preferível que
%    a decimação seja feita por último
