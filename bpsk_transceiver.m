%% Передатчик
clear all;
clc;
%% Параметры сигнала сообщения
Fm = 1e4; % Основная частота сигнала сообщения
Harm = [ 1 0.5 2 1 ]; % Частотные компоненты сигнала сообщения, представлено 4 синусоидальными рядами  
Ampl = [ 1 2 3 1 ]; % Амплитуда согласованной гармонической частотной составляющей
sampling_rate = 1/(20*max(Fm*Harm)); % Частота дискретизации системы
range = 2/min(Fm*Harm); % Диапазон выборки времени 
t = 0:sampling_rate:range; % синхронизация
%% Сигнал сообщения 
message = zeros(size(t));
for k=1:length(Harm)
    message = message + Ampl(k)*sin(2*pi*Harm(k)*Fm*t);
end
minAmpl = min(message); % Минимальная амплитуда, которая может быть использована в этой модели
%% Дискретизация
n_sample = 5;
Fs = n_sample*max(Harm*Fm); % Частота дискретизации
m_samp = zeros(size(message));
m_samp(1:1/(Fs*sampling_rate):length(t)) = message(1:1/(sampling_rate*Fs):length(t));% Сэмплированный выходной сигнал
figure(1);
plot(t,message);
grid on 
hold on
stem(t,m_samp);
xlabel('Время');
ylabel('Сигнал сообщения, сэмплированный сигнал и квантованный сигнал');
title('СИГНАЛ СООБЩЕНИЯ, СЭМПЛИРОВАННЫЙ СИГНАЛ и КВАНТОВАННЫЙ СИГНАЛ');
legend('Сигнал сообщения','Сэмплированный сигнал','Квантованный сигнал');
%% Квантование
kv_levels = 4; % Количество уровней квантования                  
quantile = (max(m_samp) - min(m_samp))/(kv_levels); % Интервал квантования
resp_levels = min(m_samp):quantile:max(m_samp); % Уровни представления
Q = zeros(size(m_samp));                                                          
for k = 1:length(resp_levels)
    values = (( m_samp>(resp_levels(k)-quantile/2) & m_samp<(resp_levels(k)+quantile/2)));   
    Q(values) = round(resp_levels(k)); % Получение квантованого сигнала выборочного сообщения
end
clear values;
stem(t,Q,'r*');
grid on
legend('Сигнал сообщения, сэмплированный сигнал и квантованный сигнал');
%% Кодирование
if min(Q) >= 0
    minAmpl = 0;                                                        
end
Q1 = Q-round(minAmpl); % Смещение отрицательных значений на положительную сторону для преобразования в двоичную выборку
Bits = de2bi(Q1(1:1/(Fs*sampling_rate):length(Q)),4,'left-msb')';          
Bits = Bits(:)'; % Генерирование двоичной последовательности дискретизированного квантованного сигнала
figure(5);
stem(Bits,'*r');
hold on;
legend('Последовательность бит в передатчике','Последовательность бит в приемнике');
%% Двоичная фазовая манипуляция (BPSK)
Fc = 1e4; % Частота несущего сигнала 
Nsamp = 10; % Образцы системы за цикл сигнала несущей 
Ncyc = 2; % Количество циклов сигнала несущей для одного битового интервала                                                         
Tbit = 0:1/(Nsamp*Fc):Ncyc/Fc; % Интервал бит
t_per = 0:1/(Nsamp*Fc):(Ncyc*length(Bits))/Fc+(length(Bits)-1)/(Nsamp*Fc); % Общее время передачи  
modulation_sig = zeros(size(t_per));                                               
l = 1;
for k = 1:length(Tbit):length(modulation_sig)
    if(Bits(l) == 1)
        modulation_sig(k:k+length(Tbit)-1) = cos(2*pi*Fc*Tbit); % Фазовая модуляция несущей для представления двоичного символа 1
    else
        modulation_sig(k:k+length(Tbit)-1) = -cos(2*pi*Fc*Tbit); % Фазовая модуляция несущей для представления двоичного символа 0
    end
    l = l+1;
end
%% AWGN - Канал
Per_sig = awgn(modulation_sig,10); % Передача модулированного сигнала несущей по каналу AWGN
figure(2);
plot(t_per,modulation_sig,'.-b',t_per,Per_sig,'r');
axis([0 3*Ncyc/Fc -2 2]);
xlabel('Время');
ylabel('Tx - передатчик сигнала и Tx - передатчик сигнала с шумом');
title('Передача сигнала и передача зашумленного сигнала');
legend('Передача сигнала','Передача зашумленного сигнала');
%% Приемник
% Фильтр
F_freq = -(Nsamp*Fc)/2:(Nsamp*Fc)/length(t_per):(Nsamp*Fc)/2-(Nsamp*Fc)/length(t_per); % Диапазон частот, используемый для визуализации БПФ сигналa 
f_pr = fft(Per_sig); % БПФ                              
figure(3);
plot(F_freq,fftshift(f_pr),F_freq,fftshift(fft(modulation_sig)),'g');
grid on;
xlabel('Частота');
ylabel('Амплитуда сигнала');
legend('Модулированный сигнал','Полученный сигнал');
title('МОДУЛИРОВАННЫЙ СИГНАЛ И ПОЛУЧЕННЫЙ СИГНАЛ');
F_rece = zeros(size(f_pr));               
FIR = (F_freq < -3*Fc | F_freq>3*Fc);
F_rece(FIR) = f_pr(FIR); % Фильтрация зашумленного сигнала в частотной области для удаления шума
F_rece(~FIR) = 0.5*f_pr(~FIR);                                                                                                                      
t_rece = ifft(F_rece); % Удаление шума из сигнала
figure(4);
plot(t_per,t_rece);
grid on;
xlabel('Время');
ylabel('Полученный сигнал');
title('ПОЛУЧЕННЫЙ СИГНАЛ ПОСЛЕ ФИЛЬТРОВАНИЯ ШУМА');
delete f_freq f_tran f_rece;
clear f_freq f_tran  f_rece;
%% Демодуляция
demod_data = zeros(size(Bits));
l = 1;
for k = 1:length(Tbit):length(t_per) % Извлечение двоичных данных из носителя с использованием метода корреляции                                            
        a = corrcoef(cos(2*pi*Fc*Tbit),t_rece(k:k+length(Tbit)-1));                                                              
        b = mean(real(a(:)));
        if b > 0.5
            demod_data(l) = 1;
        else
            demod_data(l) = 0;
        end
        l = l + 1;
end
figure(5);
stem(demod_data);
grid on;
xlabel('Позиция бит');
ylabel('Последовательность бит в приемнике и передатчике');
title('ПОСЛЕДОВАТЕЛЬНОСТЬ БИТА В ПРИЕМНИКЕ И ПЕРЕДАТЧИКЕ');
legend('Последовательность бит в приемнике','Последовательность бит в передатчике');
%% Декодирование
demod_data = reshape(demod_data,4,length(Bits)/4)';
mQ_rece = zeros(size(Q));
mQ_rece(1:1/(Fs*sampling_rate):length(Q)) = bi2de(demod_data,'left-msb')' + min(Q); % Извлечение дискретизированных квантованных данных из декодированной двоичной последовательности
%% Преобразование сигнала
F_freq = -1/(2*sampling_rate):1/(sampling_rate*length(t)):1/(2*sampling_rate)-1/(sampling_rate*length(t));                                  
F_rece = fft(mQ_rece); % БПФ полученного извлеченного квантованного сигнала с дискретизацией     
F_out = zeros(size(F_rece));                                            
figure(6);
plot(F_freq,fftshift(F_rece),F_freq,fftshift((fft(message))));
grid on;
xlabel('Частота');
ylabel('Преобразование сигнала и преобразование полученного сигнала');
title('ИСХОДНЫЙ СИГНАЛ И ПОЛУЧЕННЫЙ СИГНАЛ НА ВЫХОДЕ');
legend('Полученный сигнал','Исходный сигнал');
F_out((F_freq < -17000 | F_freq > 17000))=F_rece((F_freq < -17000 | F_freq > 17000)); % Фильтрация в частотной области для восстановления сигнала из выборки квантованных данных    
out = ifft(F_out); % Реконструированный выходной сигнал
figure(7);
plot(t,4*out,t,message,'r');
grid on;
xlabel('time');
ylabel('Сигнал сообщения и сигнал на выходе');
title('Исходный сигнал и полученный сигнал');
legend('Полученный сигнал','Исходный сигнал');
gain = 4; % Усиление мощности 
out = out*gain; % Выход после усиления мощности
figure(8);
plot(t,out);
grid on;
xlabel('Время');
ylabel('Полученный сигнал');
title('ПОЛУЧЕННЫЙ СИГНАЛ');
