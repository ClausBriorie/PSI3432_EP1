%% EXP 1 PSI3432
close all;

%% ITEM A
f = 60e9;
M=8;
c = 3e8;
lambda = c/f;
maxalpha = 0.5;
d= maxalpha*lambda;
maxchip = (M-1)*d % Comprimento máximo do chip para que não haja aliasing

%% ITEM B
theta0 = 30;
d1=lambda/4;
d2=lambda/2;
d3=3*lambda/4;
theta = -90:1:90;

B1 = 1/M*(1-exp(1i*M*2*pi*d1/lambda*(sind(theta)-sind(theta0)))) ./ ...
        (1-exp(1i*2*pi*d1/lambda*(sind(theta)-sind(theta0))));
B2 = 1/M*(1-exp(1i*M*2*pi*d2/lambda*(sind(theta)-sind(theta0)))) ./ ...
        (1-exp(1i*2*pi*d2/lambda*(sind(theta)-sind(theta0))));
B3 = 1/M*(1-exp(1i*M*2*pi*d3/lambda*(sind(theta)-sind(theta0)))) ./ ...
        (1-exp(1i*2*pi*d3/lambda*(sind(theta)-sind(theta0))));

plot(theta, abs(B1), theta, abs(B2), theta, abs(B3));shg
legend('lambda/4', 'lambda/2', '3*lambda/4');
title('Item B');

%% ITEM C
A2 = 0.5;
A1 = 1;
theta1 = 30;
theta2 = -90:1:90;

B3 = 1/M*(1-exp(1i*M*2*pi*d2/lambda*(sind(theta2)-sind(theta1)))) ./ ...
        (1-exp(1i*2*pi*d2/lambda*(sind(theta2)-sind(theta1))));
C = A1*exp(1i*2*pi*f) + B3*A2*exp(1i*2*pi*f);

figure;
plot(theta2, abs(C));shg
title('Item C');

%% ITEM D
A2 = 3*sqrt(2);
load('symbols.mat');
D2 = zeros(1,16);
theta2=60;
B3 = 1/M*(1-exp(1i*M*2*pi*d2/lambda*(sind(theta2)-sind(theta1)))) ./ ...
        (1-exp(1i*2*pi*d2/lambda*(sind(theta2)-sind(theta1))));
for k=1:1:16
    D1 = symbols(k)*exp(1i*2*pi*f) + B3*A2*exp(1i*2*pi*f);
    D2(k) = D1;
end
    
Dideal=zeros(1,16);
for k=1:1:16
    Dideal(k)=symbols(k)*exp(1i*2*pi*f);
end
figure;
scatter(real(Dideal),imag(Dideal), 'r');
hold on;
scatter(real(D2),imag(D2), 'b');
title('Item D: Diagramas de constelação ideal e para theta2 = 60°');
legend('Ideal', 'Com interferência');
grid on;
hold off;

erros = zeros(1, 181);
for theta2=-90:1:90
    D=[]; % Armazena saída para cada valor de symbols
    D3=[];
    B3 = 1/M*(1-exp(1i*M*2*pi*d2/lambda*(sind(theta2)-sind(theta1))))./ ...
            (1-exp(1i*2*pi*d2/lambda*(sind(theta2)-sind(theta1))));
    erro = 0;
    % Para cada elemento do vetor de saidas...
    for k=1:1:16
        D = symbols(k)*exp(1i*2*pi*f) + B3*A2*exp(1i*2*pi*f);
        [menorDistancia, indiceDoMinimo] = min(abs(D-Dideal));
        if Dideal(indiceDoMinimo) ~= Dideal(k)
            erro = erro + 1;
        end
    end
    erros(theta2+91) = (erro/16)*100;
end
vtheta=-90:1:90;
figure;
plot(vtheta, erros);
title('Item D - Taxa de erros em função de Theta2');
percetualDeErros = erros/16;

%ITEM E
taxaTrasm = 1e9; % em simbolos/segundo
Ta = 1e-12; % Periodo da amostragem
omega = 2*pi*f;
load('sinais.mat'); % linha <=> amostra; coluna <=> antena

sinal_saida = zeros(1, size(sinais,1));
for lin=1:size(sinais,1) % Para cada instante no vetor de sinais...
    sinal_l = 0;
    for col=1:size(sinais,2) % Para cada antena no instante...
        % Calculo o atraso a ser realizado
        tau_k = (col*d*sin(30 * ((2*pi)/360))/c);
        % Atraso este sinal
        sinal_l = sinal_l + (1/M)*sinais(lin, col)*exp(1i*omega*tau_k);
    end
    sinal_saida(lin) = sinal_l;
end
figure;
subplot(2, 1, 1)
plot(1:size(sinais, 1), sinal_saida)
title('Item E: Resultado de Beamforming Delay and Sum')
subplot(2, 1, 2)
plot(1:size(sinais, 1), sinais(:, 1))
title('Item E: Sinal recebido pela primeira antena')

%ITEM F
tam_msg = 392; % Tamanho da mensagem em bits
omega = 2*pi*f;
x_I = zeros(1, size(sinais, 1));
x_Q = zeros(1, size(sinais, 1));
for x=1:size(sinais, 1)
   x_I(x) = cos(x*omega)*sinal_saida(x);
   x_Q(x) = -sin(x*omega)*sinal_saida(x);
end
figure;
subplot(2, 1, 1)
plot(1:size(sinais, 1), x_I)
title('Item F: Produto do sinal com cosseno (x_i)')
subplot(2, 1, 2)
plot(1:size(sinais, 1), x_Q)
title('Item F: Produto do sinal com seno (x_q)')

[b, a] = butter(6, 30e9/f);
sinal_demodulado = filter(b, a, (x_I + 1i*x_Q));
figure;
plot(1:size(sinais, 1), sinal_demodulado)
title('Item F: Sinal demodulado')
