clc;
clear;
close all;

m = 8; %����� �������������� ������������������
crc = [1 0 1 0 0 1 1 0 1 0 1 1 1 1 0 0 1]; %crc16 m=8 x16+x13+x12+x11+x10+x8+x6+x5+x2+1
% crc = [1 0 1 0 1 0 1 1 1]; %crc8 m=4 x8+x7+x6+x4+x2+1
% crc = [1 1 1 0 1]; % x4+x2+x+1 m=3
% crc = [1 1 0 1]; % x3+x+1 m=4
r = length(crc)-1;
index = 1;
[H, G] = hammgen(4);
hammParam = [15, 11];
eVectors = MostSuitableEVectors(H, hammParam); % ������ �������� ���������� �������� ������ ��� ����������� ������������� �� ��������
snr_arr = -10:2:10;
snr_theor_arr = -10:10;

Pe_decode_Theor = ones(1,length(snr_theor_arr)).*(1/2^r); % ������� ������
Pb_theor = qfunc(sqrt(2*(10.^(snr_theor_arr./10)))); % ����. ����������� ������ �� ���

codewords = zeros(2^m, m+r); % ������ ������� ����
mes_3_11 = zeros(3, 11, 2^m); % ������ ������� ���� �������� �� 3 ����� �� 8 ��� � ����������� �����
hamm_3_15 = zeros(3, 15, 2^m); % ������ ���� ������������� �� ��������
modulate_3_12 = zeros(3, 12, 2^m); % ������ �������������� ������� ����
d_arr = zeros((2^m), 1); % ������ ����� ������� ����
for word = 0:2^m-1
    [~, c] = gfdeconv(de2bi(bitshift(word, r)), crc); % ����������� �����
    codewords(word+1, :) = de2bi(bitxor((bitshift(word, r)), bi2de(c)), m+r); % ������������ �������� �����
    % ������� ����� �������� �� 3 ������������������ �� 8 ��� � ����������� 3-� ����� ��� ��������� 11-�� ������ ������������������
    mes_3_11(:, :, word+1) = [reshape(codewords(word+1, :), 8, 3)' [0 0 0; 0 0 0; 0 0 0]];
    for i = 1:3
        hamm_3_15(i,:,word+1) = mod(mes_3_11(i,:,word+1)*G, 2); % ����������� �� ��������
        modulate_3_12(i,:,word+1) = hamm_3_15(i,1:12,word+1).*2-1; % �������������
    end
    d_arr(word+1) = nnz(codewords(word+1, :));
end

% ������ ������ ������ ����������� ������ ������������� CRC
d_min = min(d_arr(2:end));
P_ed = zeros(1, length(Pb_theor));
for p = Pb_theor
    tmp = 0;
    for i = d_min:(m+r)
        Ai = sum(d_arr == i);
        tmp = tmp + (Ai * p^i * (1-p)^((m+r)- i));
    end
    P_ed(index) = tmp;
    index = index + 1;
end

%������ ���� ������� ���� ��� ������� ������������� �� ��������
all_cw_for_soft = zeros(2^8, 8);
h = zeros(2^8, 12);
for i = 1:2^8
    all_cw_for_soft(i, :) = de2bi(i-1, 8);
    temp = mod([all_cw_for_soft(i, :) 0 0 0] * G, 2);
    h(i, :) = temp(1:12);
end

index = 1;
Pb = zeros(1, length(snr_arr));
Pe_decode = zeros(1, length(snr_arr));
Pe_decode_soft = zeros(1, length(snr_arr));
T = zeros(1, length(snr_arr));
tests = [100 400 1500 5000 10000 10000 10000 10000 10000 10000 10000];
for SNR = snr_arr
    % Nb - ���-�� ������� ������, Ne_decode - ���-�� ������ �������������
    % ��� ���������� ������������� �� ��������, Ne_decode_soft - ���-��
    % ������ ������������� ��� ������ ������������� �� ��������
    Nb = 0; Ne_decode = 0; Ne_decode_soft = 0; 
    SNRi = 10^(SNR/10);
    sigma = sqrt(1/(2*SNRi));
    % Nt - ���-�� ������, Nsent - ���-�� ������������ ���������
    Nt = 0; Nsent = 0;
    tic
    parfor i = 1:tests(index)
        rnd_ind = randi(2^m, 1);
        mod_mes = modulate_3_12(:, :, rnd_ind);
        a = codewords(rnd_ind, :);
        while 1
            Nt = Nt + 1;
            noisy_mes = mod_mes + sigma*randn(3, 12);
            demod_mes = noisy_mes > 0;
            hamm_synd_decode_mes = HammingSyndromeDecoder(demod_mes, H, eVectors);
            e = sum(xor(a, hamm_synd_decode_mes));
            Nb = Nb + e;
            hamm_soft_decode_mes = HammingSoftDecoder(demod_mes, all_cw_for_soft, h);
            e_soft = sum(xor(a, hamm_soft_decode_mes));
            [~, s_hammsynddecoder] = gfdeconv(hamm_synd_decode_mes, crc);
            [~, s_hammsoftdecoder] = gfdeconv(hamm_soft_decode_mes, crc);
            
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
    Pb(index) = Nb/(n*(m+r));
    Pe_decode(index) = Ne_decode/Nt;
    Pe_decode_soft(index) = Ne_decode_soft/Nt;
    T(index) = (m * Nsent)/(36 * Nt);
    index = index + 1;
end

figure(1);
semilogy(snr_theor_arr, Pe_decode_Theor, 'g', ...
         snr_theor_arr, P_ed, 'r', ...
         snr_arr, Pe_decode, 'bo', ...
         snr_arr, Pe_decode_soft, 'ko');
legend('Ped upper boound', 'Ped theor', 'Ped synd', 'Ped soft');
figure(2);
semilogy(snr_arr, T);


function eVectors = MostSuitableEVectors(H, hammParam)
% ������� ������� �������� ���������� ������� ������ ��� �������
% ��������
% ���������:
% H - ����������� ������� ���� ��������
% hammParam - ������ � ����������� ���� ��������, � ����� ������ ����� ���
%             [15, 11]
% ���������:
% eVectors - ������� �������� ���������� �������� ������ ��� �������
%            �������� � ���������� ����, �.�. ��� ����, ����� ��������
%            �������� ���������� ������ ��� �������� s = 100(2) = 4(10),
%            ����� ���������� � ����� ������� ��������� �������: 
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

function decodeCW = HammingSyndromeDecoder(demodW, H, eVectors)
% ������� ����������� ������������� ��������
% ���������:
% demodW - ���������������� ������� ����� �������� �� 3 ����� �� 12 ���,
%          �.�. ��� ������ �������� 3�12
% H - ����������� ������� ���� ��������
% eVectors - �������� ���������� ������� ������ ��� ������� ��������
% ���������:
% decodeCW - �������������� ������� ����� �������� 24 ����

    demodW = [demodW(:,:) [0 0 0; 0 0 0; 0 0 0]]; % ���������� �����, ����� ���� 15 ��� � ������ �����
    decodeCW = zeros(3, 8);
    for part = 1:3
        s = bi2de(mod(demodW(part, :)*H', 2)); % ���������� ��������
        if s > 0
            demodW(part, :) = gfadd(demodW(part, :), eVectors(s, :));
        end
        decodeCW(part, :) = demodW(part, 5:end-3);
    end
    decodeCW = reshape(decodeCW', 1, 24);
end

function decodeCW = HammingSoftDecoder(demodCW, allCW, h)
% ������� ������� ������������� �� ��������
% ���������:
% demodCW - ���������������� ������� ����� �������� �� 3 ����� �� 12 ���,
%           �.�. ��� ������ �������� 3�12
% allCW - ��� ��������� ��������� �� 8 ���
% h - ��� ��������� ������������������ � ���������� ����������� ������
%     ����� ����������� �� ��������, �� ���� ������ ������� 2^8�12
% ���������:
% decodeCW - �������������� ������� ����� �������� 24 ����

    min_d = [12 12 12];
    decodeCW = zeros(3, 8);
    for part = 1:3
        for i = 1:2^8
            min = sum(gfadd(double(demodCW(part, :)), h(i, :)));
            if min < min_d(part)
                min_d(part) = min;
                decodeCW(part, :) = allCW(i, :);
            end
        end
    end
    decodeCW = reshape(decodeCW', 1, 24);
end

