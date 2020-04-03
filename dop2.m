clc;
clear;
close all;

m = 8; %длина информационной последовательности
crc = [1 0 1 0 0 1 1 0 1 0 1 1 1 1 0 0 1]; %crc16 x16+x13+x12+x11+x10+x8+x6+x5+x2+1
r = length(crc)-1;
[H, G] = hammgen(4);
hammParam = [15, 11];
% eVectors = MostSuitableEVectors(H, hammParam); % массив наиболее подходящих векторов ошибки для синдромного декодирования по Хэммингу
eVectors = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
            0,1,0,0,0,0,0,0,0,0,0,0,0,0,0;
            0,0,0,0,1,0,0,0,0,0,0,0,0,0,0;
            0,0,1,0,0,0,0,0,0,0,0,0,0,0,0;
            0,0,0,0,0,0,0,0,1,0,0,0,0,0,0;
            0,0,0,0,0,1,0,0,0,0,0,0,0,0,0;
            0,0,0,0,0,0,0,0,0,0,1,0,0,0,0;
            0,0,0,1,0,0,0,0,0,0,0,0,0,0,0;
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,1;
            0,0,0,0,0,0,0,0,0,1,0,0,0,0,0;
            0,0,0,0,0,0,0,1,0,0,0,0,0,0,0;
            0,0,0,0,0,0,1,0,0,0,0,0,0,0,0;
            0,0,0,0,0,0,0,0,0,0,0,0,0,1,0;
            0,0,0,0,0,0,0,0,0,0,0,1,0,0,0;
            0,0,0,0,0,0,0,0,0,0,0,0,1,0,0];
snr_arr = -10:2:10;
snr_theor_arr = -10:10;

Pe_decode_Theor = ones(1,length(snr_theor_arr)).*(1/2^r); % верхняя оценка
Pb_theor = qfunc(sqrt(2*(10.^(snr_theor_arr./10)))); % теор. вероятность ошибки на бит

codewords = zeros(2^m, m+r); % массив кодовых слов
mes_3_11 = zeros(3, 11, 2^m); % массив кодовых слов разбитых на 3 части по 8 бит с добавлением нулей
hamm_3_15 = zeros(3, 15, 2^m); % массив слов закодированых по Хэммингу
modulate_3_12 = zeros(3, 12, 2^m); % массив модулированных кодовых слов
d_arr = zeros((2^m), 1); % массив весов кодовых слов
for word = 0:2^m-1
    [~, c] = gfdeconv(de2bi(bitshift(word, r)), crc); % контрольная сумма
    codewords(word+1, :) = de2bi(bitxor((bitshift(word, r)), bi2de(c)), m+r); % формирование кодового слова
    % кодовые слова разбитые на 3 последовательности по 8 бит с добавлением 3-х нулей для получения 11-ти битной последовательности
    mes_3_11(:, :, word+1) = [reshape(codewords(word+1, :), 8, 3)' [0 0 0; 0 0 0; 0 0 0]];
    for i = 1:3
        hamm_3_15(i,:,word+1) = mod(mes_3_11(i,:,word+1)*G, 2); % кодирование по Хэммингу
        modulate_3_12(i,:,word+1) = hamm_3_15(i,1:12,word+1).*2-1; % модулирование
    end
    d_arr(word+1) = nnz(codewords(word+1, :));
end

% Расчет точной оценки вероятности ошибки декодирования CRC
d_min = min(d_arr(2:end));
P_ed = zeros(1, length(Pb_theor));
ind = 1;
for p = Pb_theor
    tmp = 0;
    for i = d_min:(m+r)
        Ai = sum(d_arr == i);
        tmp = tmp + (Ai * p^i * (1-p)^((m+r)- i));
    end
    P_ed(ind) = tmp;
    ind = ind + 1;
end

%Расчет всех кодовых слов для мягкого декодирования по Хэммингу
all_cw_for_soft = zeros(2^8, 8);
h = zeros(2^8, 12);
for i = 1:2^8
    all_cw_for_soft(i, :) = de2bi(i-1, 8);
    temp = mod([all_cw_for_soft(i, :) 0 0 0] * G, 2);
    h(i, :) = temp(1:12)*2-1;
end

ind = 1;
Pb = zeros(1, length(snr_arr));
Pb_synd = zeros(1, length(snr_arr));
Pb_soft = zeros(1, length(snr_arr));
Pe_decode = zeros(1, length(snr_arr));
Pe_decode_soft = zeros(1, length(snr_arr));
T = zeros(1, length(snr_arr));
tests = [100 400 1500 5000 10000 10000 10000 10000 10000 10000 10000];
for SNR = snr_arr
    % Nb - кол-во битовых ошибок, Ne_decode - кол-во ошибок декодирования
    % при синдромном декодировании по Хэммингу, Ne_decode_soft - кол-во
    % ошибок декодирования при мягком декодировании по Хэммингу
    Nb = 0; Nb_synd = 0; Nb_soft = 0; Ne_decode = 0; Ne_decode_soft = 0; 
    SNRi = 10^(SNR/10);
    sigma = sqrt(1/(2*SNRi));
    % Nt - кол-во тестов, Nsent - кол-во отправленных сообщений
    Nt = 0; Nsent = 0;
    tic
    for i = 1:tests(ind)
        rnd_ind = randi(2^m, 1);
        mod_mes = modulate_3_12(:, :, rnd_ind);
        a = codewords(rnd_ind, :);
        while 1
            Nt = Nt + 1;
            noisy_mes = mod_mes + sigma*randn(3, 12);
            demod_mes = noisy_mes > 0;
            [hamm_synd_decode_mes, nErrBits, v_synd] = HammingSyndromeDecoder(noisy_mes, H, eVectors, mod_mes);
            Nb_synd = Nb_synd + v_synd;
            Nb = Nb + nErrBits;
            [hamm_soft_decode_mes, v_soft] = HammingSoftDecoder(noisy_mes, all_cw_for_soft, h, mod_mes);
            Nb_soft = Nb_soft + v_soft;
            [~, s_hammsynddecoder] = gfdeconv(hamm_synd_decode_mes, crc);
            e = sum(xor(a, hamm_synd_decode_mes));
            [~, s_hammsoftdecoder] = gfdeconv(hamm_soft_decode_mes, crc);
            e_soft = sum(xor(a, hamm_soft_decode_mes));

            
            if (bi2de(s_hammsoftdecoder) == 0) && (e_soft > 0)
                Ne_decode_soft = Ne_decode_soft + 1;
            end
            
            if bi2de(s_hammsynddecoder) == 0
                if e > 0
                    Ne_decode = Ne_decode + 1;
                end
                Nsent = Nsent + 1;
                break;
            end
        end
        
    end
    toc
    Pb(ind) = Nb/(Nt * 36);
    Pb_synd(ind) = Nb_synd/(Nt * 36);
    Pb_soft(ind) = Nb_soft/(Nt * 36);
    Pe_decode(ind) = Ne_decode/Nt;
    Pe_decode_soft(ind) = Ne_decode_soft/Nt;
    T(ind) = (m * Nsent)/(36 * Nt);
    ind = ind + 1;
end

figure(1);
semilogy(snr_theor_arr, Pe_decode_Theor, 'g', ...
         snr_theor_arr, P_ed, 'r', ...
         snr_arr, Pe_decode, 'bo', ...
         snr_arr, Pe_decode_soft, 'ko');
legend('Ped upper boound', 'Ped theor', 'Ped synd', 'Ped soft');
figure(2);
semilogy(snr_theor_arr, Pb_theor, ...
         snr_arr, Pb, 'o', ...
         snr_arr, Pb_synd, '*', ...
         snr_arr, Pb_soft, 'x')
legend('Pb_theor', 'Pb', 'Pb_synd', 'Pb_soft')
figure(3);
semilogy(snr_arr, T);


function eVectors = MostSuitableEVectors(H, hammParam)
% Функция находит наиболее подходящие вектора ошибок для каждого
% синдрома
% Аргументы:
% H - проверочная матрица кода Хэмминга
% hammParam - массив с параметрами кода Хэмминга, в нашем случае имеет вид
%             [15, 11]
% Результат:
% eVectors - таблица наиболее подходящих векторов ошибки для каждого
%            синдрома в десятичном виде, т.е. для того, чтобы получить
%            наиболее подходящий вектор для синдрома s = 100(2) = 4(10),
%            нужно обратиться к этому массиву слудующим образом: 
%            eVectors(4, :)

    s = 1:bi2de(ones(1, hammParam(1)-hammParam(2)));
    eVectors = zeros(length(s), hammParam(1));
    for i = s
        tmp = [];
        for e = 1:bi2de(ones(1, hammParam(1)))
            if i == bi2de(mod(de2bi(e, hammParam(1))*H', 2))
                tmp = [tmp e];
            end
        end
        min = sum(ones(1, hammParam(1)));
        for j = 1:length(tmp)
            if sum(de2bi(tmp(j))) < min
                min = sum(de2bi(tmp(j)));
                eVectors(i, :) = de2bi(tmp(j), hammParam(1));
            end
        end
    end
end

function [decodeCW, nErrBits, v] = HammingSyndromeDecoder(noisy_mes, H, eVectors, modCW)
% Функция синдромного декодирования Хэмминга
% Аргументы:
% noisy_mes - зашумленное сообщение разбитое на 3 части по 12 бит,
%             т.е. это массив размером 3х12
% H - проверочная матрица кода Хэмминга
% eVectors - наиболее подходящие вектора ошибок для каждого синдрома
% modCW - модулированное кодовое слово
% Результаты:
% decodeCW - декодированное кодовое слово размером 24 бита
% nErrBits - количество битовых ошибок
% v - количество битовых ошибок после кода Хэмминга
    
    nErrBits = sum(sum(xor(noisy_mes > 0, modCW > 0)));
    demodW = [noisy_mes(:,:)>0 [0 0 0; 0 0 0; 0 0 0]]; % добавление нулей, чтобы было 15 бит в каждой части
    decode = zeros(3, 12);
    for part = 1:3
        s = bi2de(mod(demodW(part, :)*H', 2)); % вычисление синдрома
        if s > 0
            demodW(part, :) = gfadd(demodW(part, :), eVectors(s, :));
        end
        decode(part, :) = demodW(part, 1:end-3);
    end
    v = sum(sum(xor(decode, modCW > 0)));
    decodeCW = reshape(decode(:, 5:end)', 1, 24);
end

function [decodeCW, v] = HammingSoftDecoder(noisy_mes, allCW, h, modCW)
% Функция мягкого декодирования по Хэммингу
% Аргументы:
% noisy_mes - зашумленное сообщение разбитое на 3 части по 12 бит,
%             т.е. это массив размером 3х12
% allCW - все возможные сообщения по 8 бит
% h - все возможные последовательности с удаленными добавочными нулями
%     после кодирования по Хэммингу и модулирования, т.е. массив размера 2^8х12
% modCW - модулированное кодовое слово
% Результаты:
% decodeCW - декодированное кодовое слово размером 24 бита
% v - количество битовых ошибок после кода Хэмминга

    min_d = [100 100 100];
    decodeCW = zeros(3, 8);
    cor = zeros(3, 12);
    for part = 1:3
        for i = 1:2^8
            min = sqrt(sum((noisy_mes(part, :) - h(i, :)).^2));
            if min < min_d(part)
                min_d(part) = min;
                decodeCW(part, :) = allCW(i, :);
                cor(part, :) = h(i, :);
            end
        end
    end
    v = sum(sum(xor(cor > 0, modCW > 0)));
    decodeCW = reshape(decodeCW', 1, 24);
end
