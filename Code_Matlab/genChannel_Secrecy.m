function [channels, sys] = genChannel_Secrecy(sys)
    % Generate channels for SIM-ISAC Secrecy Paper
    Mx = sys.Mx; Mz = sys.Mz; M = sys.M; Nt = sys.Nt; K = sys.K;
    lambda = sys.lambda; df = sys.df; L = sys.L; d_layer = sys.d_layer;
    C0 = sys.C0; n_bar = sys.n_bar;

    %% 1. Transmission Matrices (Rayleigh-Sommerfeld)
    Psi = zeros(M, M);
    for m = 1:M
        mx = mod(m-1, Mx); mz = floor((m-1)/Mx);
        for mp = 1:M
            mpx = mod(mp-1, Mx); mpz = floor((mp-1)/Mx);
            dist_sq = (mx-mpx)^2 * df^2 + (mz-mpz)^2 * df^2;
            d_mm_p = sqrt(d_layer^2 + dist_sq);
            cos_theta = d_layer / d_mm_p;
            Psi(m, mp) = (sys.dx * sys.dy * cos_theta / d_mm_p) * ...
                         (1/(2*pi*d_mm_p) - 1i/lambda) * exp(1i*2*pi*d_mm_p/lambda);
        end
    end
    channels.Psi = Psi;

    % BS to 1st Layer (H_BS_SIM)
    H_BS_SIM = zeros(M, Nt);
    for m = 1:M
        mx = mod(m-1, Mx); mz = floor((m-1)/Mx);
        for n = 1:Nt
            nx = n - (Nt+1)/2;
            dist_sq = (mx*df - nx*df)^2 + (mz*df)^2;
            d_mn = sqrt(d_layer^2 + dist_sq);
            cos_theta = d_layer / d_mn;
            H_BS_SIM(m, n) = (sys.dx * sys.dy * cos_theta / d_mn) * ...
                             (1/(2*pi*d_mn) - 1i/lambda) * exp(1i*2*pi*d_mn/lambda);
        end
    end
    channels.H = H_BS_SIM;

    %% 2. Spatial Correlation Matrix R (M x M)
    R = zeros(M, M);
    for m = 1:M
        mx = mod(m-1, Mx); mz = floor((m-1)/Mx);
        for mp = 1:M
            mpx = mod(mp-1, Mx); mpz = floor((mp-1)/Mx);
            d_m_mp = sqrt((mx-mpx)^2 + (mz-mpz)^2) * df;
            if d_m_mp == 0
                R(m, mp) = 1;
            else
                R(m, mp) = sinc(2 * d_m_mp / lambda);
            end
        end
    end
    R = R + 1e-9 * eye(M);
    R_sqrt = chol(R, 'lower');
    channels.R = R;
    channels.R_sqrt = R_sqrt;   % Lưu thêm để tiện dùng

    %% 3. SIM to UEs and Eves Channels - VỊ TRÍ UE NGẪU NHIÊN
    T_SIM = sys.T_SIM;
    
    % Tham số random placement
    D_center = 80;      % Khoảng cách trung bình từ SIM (m)
    D_spread = 30;      % Độ dao động khoảng cách (±15m)

    channels.h = zeros(M, K);
    channels.v = zeros(M, K);
    
    % Lưu thêm thông tin vị trí (rất hữu ích sau này)
    channels.x_UE = zeros(K,1);
    channels.y_UE = zeros(K,1);
    channels.z_UE = zeros(K,1);
    channels.d_ue  = zeros(K,1);

    for k = 1:K
        % === Sinh vị trí UE ngẫu nhiên ===
        phi_k   = deg2rad(-60) + rand() * deg2rad(120);     % Góc ngẫu nhiên -60° đến +60°
        dist_k  = D_center + (rand()-0.5) * D_spread;       % Khoảng cách ngẫu nhiên
        
        channels.x_UE(k) = dist_k * cos(phi_k);
        channels.y_UE(k) = dist_k * sin(phi_k);
        channels.z_UE(k) = 1.65;                            % Chiều cao UE tiêu chuẩn
        
        % Tính khoảng cách từ SIM đến UE
        dx = channels.x_UE(k) - 0;
        dy = channels.y_UE(k) - 0;
        dz = channels.z_UE(k) - T_SIM;
        
        d_ue = sqrt(dx^2 + dy^2 + dz^2);
        channels.d_ue(k) = d_ue;
        
        alpha_ue = C0 * (d_ue ^ (-n_bar));

        % === Sinh kênh h (UE) ===
        channels.h(:, k) = sqrt(alpha_ue) * R_sqrt * (randn(M, 1) + 1i*randn(M, 1)) / sqrt(2);

        % === Eve (cao hơn UE 2m) ===
        dz_eve = (channels.z_UE(k) + 2) - T_SIM;
        d_eve = sqrt(dx^2 + dy^2 + dz_eve^2);
        alpha_eve = C0 * (d_eve ^ (-n_bar));
        
        channels.v(:, k) = sqrt(alpha_eve) * R_sqrt * (randn(M, 1) + 1i*randn(M, 1)) / sqrt(2);
    end

    %% 4. Hướng điều khiển cho mục tiêu cảm biến
    channels.a = steering_vector_Secrecy(sys, sys.theta_c, sys.phi_c);
end

%% ====================== SUB FUNCTION ======================
function a = steering_vector_Secrecy(sys, theta, phi)
    Mx = sys.Mx; Mz = sys.Mz; M = sys.M;
    ax = exp(-1i * pi * sin(theta) * sin(phi) * (0:Mx-1).');
    az = exp(-1i * pi * cos(theta) * (0:Mz-1).');
    a = kron(ax, az) / sqrt(M);
end