% SIM_ISAC_config_Secrecy.m
clear

%% Tham số mô phỏng
sys.n_monte     = 50;     % Số lượt mô phỏng Monte Carlo
sys.itermax     = 30;     % Số vòng lặp tối đa
sys.itermax_init= 20;
sys.tol         = 1e-6;  % Ngưỡng hội tụ
sys.ga_iter     = 1;      % Số bước nhảy GA
sys.mu_init     = 0.5;    % Bước nhảy ban đầu

%% Thông số hệ thống
sys.K   = 4;              % Số người dùng (và Eves)
sys.Nt  = 32;             % Số anten BS
sys.fc  = 28e9;           % Tần số sóng mang
c       = 3e8;            % Vận tốc ánh sáng
sys.lambda = c / sys.fc;  % Bước sóng

    % SIM param
    sys.M       = 49;                     % Tổng số phần tử mỗi lớp
    sys.Mx      = 7;                      % Số phần tử theo trục X
    sys.Mz      = 7;                      % Số phần tử theo trục Z
    sys.L       = 3;                      % Số lớp SIM
    sys.dx      = sys.lambda / 2;
    sys.dy      = sys.lambda / 2;
    sys.df      = sys.lambda / 2;         % Khoảng cách giữa các meta-atom
    sys.T_total = 5 * sys.lambda;         % Tổng độ dày SIM
    sys.d_layer = sys.T_total / sys.L;    % Khoảng cách giữa các lớp
    
    %Thông số hình học
    sys.H_BS        = 10;                 % Độ cao BS (m)
    sys.T_SIM       = 5;                  % Độ cao SIM (m)
    sys.d_UE_spacing= 2;                  % Khoảng cách giữa các UE (m)

    % Channel param
    sys.C0_dB = -32;        % Độ suy hao tại khoảng cách tham chiếu 1 mét
    sys.C0 = 10^(sys.C0_dB/10);         
    sys.n_bar = 3.5;        % Hệ số suy hao môi trường
    sys.Rmin = 0.5;         % Tốc độ truyền dẫn tối thiểu (QoS) cho mỗi người dùng phải đạt


% Power Setup
sys.Pmax_dBm = 30;                    % Công suất phát tối đa BS
sys.Pmax     = 10^(sys.Pmax_dBm/10);
sys.rho      = 0.35;                  % Hiệu suất PA
sys.P_meta   = 1.5e-3;                % Công suất mỗi meta-atom (W)
sys.P_ctrl   = 1;                     % Công suất controller SIM (W)
sys.P0       = 4.5;

%% Sensing Parameters
sys.Gamma_dBm = -10;                  % Ngưỡng Sensing
sys.Gamma     = 10^(sys.Gamma_dBm/10);
sys.theta_c   = deg2rad(30);          % Góc ngẩng mục tiêu
sys.phi_c     = deg2rad(45);          % Góc phương vị mục tiêu

%% Noise Parameters
sys.sigma2_u_dBm = -90;               % Nhiễu tại UE
sys.sigma2_e_dBm = -90;               % Nhiễu tại Eve
sys.sigma2_u     = 10^(sys.sigma2_u_dBm/10);
sys.sigma2_e     = 10^(sys.sigma2_e_dBm/10);

fprintf('--- Secrecy SIM-ISAC System Configured ---\n');
