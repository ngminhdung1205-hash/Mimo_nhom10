% Khởi tạo các tham số
M = 4; alphaM = 2; betaM = 1;
d = 1; alpha = 2;
Mt = 2; Mr = 2;

PN0dB = 0:2:25; PN0 = 10.^(PN0dB/10);
N0 = 0.1; P = PN0*N0;
Ntry = 1*10^3; Ns = 1*10^3;

%------- 1: Kiểm chứng khôi phục tín hiệu (H không đổi, nhiễu = 0)---------
H_const = [0.7 0.3; 0.4 0.6];
X_test  = [1+1i, 1-1i];

% a. Kiểm chứng bằng phương pháp STBC (Alamouti 2x2)
r_test = zeros(1,2); z_test = zeros(1,2); s_de_STBC_test = zeros(1,2);

% Đầu thu (nhiễu = 0)
r_test(1) =  H_const(1,1)*X_test(1) + H_const(1,2)*X_test(2);
r_test(2) = -H_const(1,1)*conj(X_test(2)) + H_const(1,2)*conj(X_test(1));
z_test(1) =  H_const(2,1)*X_test(1) + H_const(2,2)*X_test(2);
z_test(2) = -H_const(2,1)*conj(X_test(2)) + H_const(2,2)*conj(X_test(1));

% Giải điều chế / Kết hợp Alamouti
temp_test = abs(H_const(1,1))^2 + abs(H_const(1,2))^2 + abs(H_const(2,1))^2 + abs(H_const(2,2))^2;
s_de_STBC_test(1) = (conj(H_const(1,1))*r_test(1) + H_const(1,2)*conj(r_test(2)) + conj(H_const(2,1))*z_test(1) + H_const(2,2)*conj(z_test(2)))/temp_test;
s_de_STBC_test(2) = (conj(H_const(1,2))*r_test(1) - H_const(1,1)*conj(r_test(2)) + conj(H_const(2,2))*z_test(1) - H_const(2,1)*conj(z_test(2)))/temp_test;

% b. Kiểm chứng bằng phương pháp SDM (SVD Precoding & Combining, slide 11 & 14)
[U_sdm, S_sdm, V_sdm] = svd(H_const); % Phân tích SVD: H = U*S*V'
X_sdm_in = X_test.';                  % Vector cột tín hiệu phát (dùng .' để không lấy liên hợp)
X_pre    = V_sdm * X_sdm_in;          % Tiền mã hóa tại bộ phát (x = V*x_tilde)
Y_sdm    = H_const * X_pre;           % Tín hiệu qua kênh H (nhiễu = 0)
Y_tilde  = U_sdm' * Y_sdm;            % Kết hợp tại bộ thu
s_de_SDM_test = (S_sdm \ Y_tilde).';  % Chia cho độ lợi kênh con sigma_i để khôi phục X

disp('1. Kết quả khi nhiễu cộng bằng 0)');
disp('Tín hiệu phát ban đầu X:'); disp(X_test);
disp('Tín hiệu thu STBC (r, z):'); disp(r_test); disp(z_test);
disp('Tín hiệu khôi phục STBC:'); disp(s_de_STBC_test);
disp('Tín hiệu tiền mã hóa SDM (X_pre) và tại đầu thu (Y_sdm):'); disp(X_pre.'); disp(Y_sdm.');
disp('Tín hiệu khôi phục SDM (SVD):'); disp(s_de_SDM_test);

%--------- 2. Khảo sát hiệu năng BER (so sánh nhiều loại phân tập)---------
BER1_1_simul = zeros(1, length(PN0dB)); % Sửa lỗi zeros(N) -> zeros(1, N)
BER2_1_simul = zeros(1, length(PN0dB)); % Bổ sung mô phỏng 2x1 để so sánh phân tập
BER2_2_simul = zeros(1, length(PN0dB));

for i = 1:length(PN0dB)
    Num_Err_1x1 = 0;
    Num_Err_2x1 = 0;
    Num_Err     = 0;
    for j = 1:Ntry
        %--------------------- Đầu phát -----------------------
        [s,Bit_data,Point_con] = a__Modulation_QPSK(Ns);

        %------------------ Kênh truyền fading ----------------
        h11 = sqrt(1/2) * (randn(1,1) + 1i*randn(1,1));
        h12 = sqrt(1/2) * (randn(1,1) + 1i*randn(1,1));
        h21 = sqrt(1/2) * (randn(1,1) + 1i*randn(1,1));
        h22 = sqrt(1/2) * (randn(1,1) + 1i*randn(1,1));

        r_1x1 = zeros(1,Ns); s_de_1x1 = zeros(1,Ns);
        r     = zeros(1,Ns); s_de_2x1 = zeros(1,Ns);
        z     = zeros(1,Ns); s_de     = zeros(1,Ns);

        for k = 1:2:Ns
            %-------------------- AWGN ------------------------
            % Đưa vào trong vòng k để mỗi cặp ký hiệu nhận nhiễu độc lập
            n1 = sqrt(N0/2) * (randn(1,1) + 1i*randn(1,1));
            n2 = sqrt(N0/2) * (randn(1,1) + 1i*randn(1,1));
            n3 = sqrt(N0/2) * (randn(1,1) + 1i*randn(1,1));
            n4 = sqrt(N0/2) * (randn(1,1) + 1i*randn(1,1));

            %---------Tín hiệu thu & tách sóng 1x1 (không phân tập) -------
            r_1x1(k)   = sqrt(P(i)*d^(-alpha)) * h11 * s(k)   + n1;
            r_1x1(k+1) = sqrt(P(i)*d^(-alpha)) * h11 * s(k+1) + n2;
            temp_1x1   = sqrt(P(i)*d^(-alpha)) * (abs(h11)^2);
            s_de_1x1(k)   = (conj(h11)*r_1x1(k))   / temp_1x1;
            s_de_1x1(k+1) = (conj(h11)*r_1x1(k+1)) / temp_1x1;

            %-------- Tín hiệu thu 2x1 & 2x2 (Alamouti STBC) --------------
            r(k)   =  sqrt(P(i)*d^(-alpha)) * h11 * s(k)         + sqrt(P(i)*d^(-alpha)) * h12 * s(k+1)       + n1;
            r(k+1) = -sqrt(P(i)*d^(-alpha)) * h11 * conj(s(k+1)) + sqrt(P(i)*d^(-alpha)) * h12 * conj(s(k))   + n2;
            z(k)   =  sqrt(P(i)*d^(-alpha)) * s(k)*h21           + sqrt(P(i)*d^(-alpha)) * s(k+1)*h22         + n3;
            z(k+1) = -sqrt(P(i)*d^(-alpha)) * conj(s(k+1))*h21   + sqrt(P(i)*d^(-alpha)) * conj(s(k))*h22     + n4;

            %------- Tín hiệu giải điều chế 2x1 (phân tập phát) -----------
            temp_2x1 = sqrt(P(i)*d^(-alpha)) * (abs(h11)^2 + abs(h12)^2);
            s_de_2x1(k)   = (conj(h11)*r(k) + h12*conj(r(k+1)))/temp_2x1;
            s_de_2x1(k+1) = (conj(h12)*r(k) - h11*conj(r(k+1)))/temp_2x1;

            %-------- Tín hiệu giải điều chế 2x2 (phân tập phát-thu) ------
            temp = sqrt(P(i)*d^(-alpha)) * (abs(h11)^2 + abs(h12)^2 + abs(h21)^2 + abs(h22)^2);
            s_de(k)   = (conj(h11)*r(k) + h12*conj(r(k+1)) + conj(h21)*z(k) + h22*conj(z(k+1)))/temp;
            s_de(k+1) = (conj(h12)*r(k) - h11*conj(r(k+1)) + conj(h22)*z(k) - h21*conj(z(k+1)))/temp;
        end

        %---------------------- Giải mã ----------------------
        Bit_de_1x1 = a__Demodulation_QPSK(s_de_1x1, Point_con);
        Bit_de_2x1 = a__Demodulation_QPSK(s_de_2x1, Point_con);
        Bit_de     = a__Demodulation_QPSK(s_de, Point_con);

        %---------------------- Error ------------------------
        [Num_1x1, ~] = biterr(Bit_data, Bit_de_1x1);
        [Num_2x1, ~] = biterr(Bit_data, Bit_de_2x1);
        [Num, ~]     = biterr(Bit_data, Bit_de);

        Num_Err_1x1 = Num_Err_1x1 + Num_1x1;
        Num_Err_2x1 = Num_Err_2x1 + Num_2x1;
        Num_Err     = Num_Err     + Num;
    end
    BER1_1_simul(i) = Num_Err_1x1/(length(Bit_data)*Ntry);
    BER2_1_simul(i) = Num_Err_2x1/(length(Bit_data)*Ntry);
    BER2_2_simul(i) = Num_Err/(length(Bit_data)*Ntry);
end

%---------------------- Lý thuyết BER ---------------------------
phiM = PN0 * betaM * d^(-alpha);
BER1_asy = (1/(2*log2(M)))*alphaM*betaM^-1*d^(1*alpha)*PN0.^-1;
BER2_asy = (3/(4*log2(M)))*alphaM*betaM^-2*d^(2*alpha)*PN0.^-2;
BER1 = alphaM/(2*log2(M))*(1-sqrt(0.5*phiM ./ (1 + 0.5*phiM)));
BER2 = alphaM/(2*log2(M)) * (1 - sqrt(phiM./(2+phiM)) - 1./(2+phiM) .* sqrt(phiM./(2+phiM)));
BER2_2 = zeros(1,length(PN0dB));
for i = 1:length(PN0dB)
    Sum = 0;
    for k = 0:Mt*Mr-1
        temp = gamma(k+1/2)/factorial(k)*(2/(2+phiM(i)))^k*sqrt(phiM(i)/(2+phiM(i)));
        Sum = Sum + temp;
    end
    BER2_2(i) = alphaM/(2*log2(M))*(1-Sum/sqrt(pi));
end

% --------3. Khảo sát tốc độ đạt được (MIMO SDM vs STBC vs SISO)-----------
C_SISO    = zeros(1, length(PN0dB));
C_STBC    = zeros(1, length(PN0dB));
C_SDM_EPA = zeros(1, length(PN0dB));
C_SDM_WF  = zeros(1, length(PN0dB));

sigma_h = svd(H_const); % Các trị riêng của kênh H_const

for i = 1:length(PN0dB)
    P_tot = P(i) * d^(-alpha); % Tổng công suất phát có tính suy hao đường truyền

    % 1. Dung lượng SISO 1x1 (lấy nhánh h11 của H_const)
    C_SISO(i) = log2(1 + P_tot * abs(H_const(1,1))^2 / N0);

    % 2. Dung lượng MIMO 2x2 STBC (bậc ghép kênh = 1)
    C_STBC(i) = log2(1 + (P_tot / Mt) * sum(abs(H_const(:)).^2) / N0);

    % 3. Dung lượng MIMO 2x2 SDM phân bổ công suất đều
    C_SDM_EPA(i) = sum(log2(1 + (P_tot / Mt) * (sigma_h.^2) / N0));

    % 4. Dung lượng MIMO 2x2 SDM tối ưu bằng Waterfilling (slide 17 & 18)
    inv_g = N0 ./ (sigma_h.^2);
    inv_g_sorted = sort(inv_g, 'ascend');
    for m = length(sigma_h):-1:1
        mu = (P_tot + sum(inv_g_sorted(1:m))) / m;
        if mu > inv_g_sorted(m)
            break;
        end
    end
    q_wf = max(0, mu - inv_g);
    C_SDM_WF(i) = sum(log2(1 + q_wf .* (sigma_h.^2) / N0));
end

%------------------------Vẽ đồ thị BER và Dung lượng------------------------------
figure(1); clf;

% Hình bên trái: BER
subplot(1,2,1);
semilogy(PN0dB,BER1,'b','linewidth',1.2); hold on;
semilogy(PN0dB,BER1_1_simul,'bo','linewidth',1.2,'markersize',5);
semilogy(PN0dB,BER2,'r','linewidth',1.2);
semilogy(PN0dB,BER2_1_simul,'r^','linewidth',1.2,'markersize',5);
semilogy(PN0dB,BER2_2,'g','linewidth',1.2);
semilogy(PN0dB,BER2_2_simul,'ko','linewidth',1.2,'markersize',5);
hold off;
xlabel('SNR (dB)'); ylabel('BER'); grid on;
title('Hiệu năng BER các cấu hình phân tập MIMO');
legend('Exact 1x1', 'Simul 1x1', 'Exact 2x1', 'Simul 2x1', 'Exact 2x2', 'Simulation 2x2', 'Location', 'SouthWest');

% Hình bên phải: Dung lượng
subplot(1,2,2);
plot(PN0dB, C_SISO, 'b-s', 'linewidth', 1.2); hold on;
plot(PN0dB, C_STBC, 'r-^', 'linewidth', 1.2);
plot(PN0dB, C_SDM_EPA, 'g-o', 'linewidth', 1.2);
plot(PN0dB, C_SDM_WF, 'k--', 'linewidth', 1.5);
hold off;
xlabel('SNR (dB)'); ylabel('Dung lượng kênh (bit/s/Hz)'); grid on;
title('So sánh dung lượng với kênh H cố định');
legend('SISO 1x1', 'MIMO 2x2 STBC', 'MIMO 2x2 SDM (Equal Power)', 'MIMO 2x2 SDM (Waterfilling)', 'Location', 'NorthWest');

% Các hàm điều chế và giải điều chế
function [si,Bit_data_NU,Point_con] = a__Modulation_QPSK(Ns1)

% QPSK Modulation
aIQ = 1/sqrt(2);
Point_con = [aIQ+1i*aIQ  aIQ-1i*aIQ  -aIQ+1i*aIQ  -aIQ-1i*aIQ];
Bit_data_NU = randsrc(1,Ns1*2,[0,1]);
qamdata = bi2de(reshape(Bit_data_NU,Ns1,2),'left-msb')';
si = genqammod(qamdata, Point_con);

end

function [data] = a__Demodulation_QPSK(Zn, Point_con)
ri = genqamdemod(Zn,Point_con);
data = de2bi(ri',2,'left-msb');
data = [data(:,1)', data(:,2)'];

end