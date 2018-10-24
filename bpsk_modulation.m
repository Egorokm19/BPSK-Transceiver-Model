clear all;
close all;

f = 2; % частота синусоидального сигнала 
fs = 10000; % периода выборки синусоидального сигнала
t = 0:1/fs:1; % время разбиения на сегменты 1 / fs
%% Установка фазовых сдвигов для различных сигналов BPSK
p1 = 0;
p2 = pi;
% Получение числовых значений битов для модуляции
N = 8;
% Генерируем случайный сигнал
bit_stream=round(rand(1,N));
% Распределение динамических переменных
time = [];
digital_signal = [];
BPSK = [];
carrier_signal = [];

%% Генерация сигнала
for ii = 1:1:N
% исходный цифровой сигнал
if bit_stream(ii) == 0
bit = zeros(1,length(t));
else
bit = ones(1,length(t));
end
digital_signal = [digital_signal bit];

% Генерирование сигнала BPSK
if bit_stream(ii) == 0
bit = sin(2*pi*f*t+p1);
else
bit = sin(2*pi*f*t+p2);
end
BPSK = [BPSK bit];

% Генерация несущей 
carrier = sin(2*f*t*pi);
carrier_signal = [carrier_signal carrier];

time = [time t];
t = t + 1;

end

subplot(3,1,1);
plot(time,digital_signal,'r');
grid on;
title('Цифровой сигнал')
axis([0 time(end) -0.5 1.5]);

subplot(3,1,2);
plot(time,BPSK);
title('BPSK')
grid on;
axis tight;

subplot(3,1,3);
plot(time,carrier_signal);
title('Сигнал несущей')
grid on;
axis tight;
