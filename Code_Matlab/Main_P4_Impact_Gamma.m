% =========================================================================
% File: Main_P4_Impact_Gamma.m
% Mục đích: Khảo sát ảnh hưởng của Ngưỡng Sensing (Sự đánh đổi ISAC)
%           So sánh Proposed (tối ưu W, Phi) và Baseline (Random Phase)
% =========================================================================
clear; clc; close all;
%% 1. KHỞI TẠO CẤU HÌNH GỐC
SIM_ISAC_config_Secrecy; 
sys_base = sys;
fprintf('\n>>> RUNNING PART 4: IMPACT OF SENSING THRESHOLD (%d Monte Carlo) <<<\n', sys.n_monte);

%% 2. VÒNG LẶP MONTE CARLO
Gamma_dBm_list = -20:5:10; 
WSR_vs_Gamma_opt_accum = zeros(length(Gamma_dBm_list), 1);
WSR_vs_Gamma_rand_accum = zeros(length(Gamma_dBm_list), 1); % Thêm biến Baseline

for m = 1:sys_base.n_monte
    fprintf('Part 4 - Monte Carlo [%d / %d]\n', m, sys_base.n_monte);
    [channels, ~] = genChannel_Secrecy(sys_base); 
    
    % Khởi tạo pha SIM ngẫu nhiên dùng chung
    Phi_init = cell(sys_base.L, 1);
    for l = 1:sys_base.L, Phi_init{l} = diag(exp(1i * 2 * pi * rand(sys_base.M, 1))); end
    
    for idx = 1:length(Gamma_dBm_list)
        sys = sys_base;
        sys.Gamma_dBm = Gamma_dBm_list(idx);
        sys.Gamma = 10^(sys.Gamma_dBm / 10); 
        
        % Khởi tạo W theo MRT (dùng chung)
        Theta_init = eye(sys.M);
        for l = 1:sys_base.L, Theta_init = Phi_init{l} * Theta_init; if l < sys.L, Theta_init = channels.Psi * Theta_init; end, end
        H_eff_init = Theta_init * channels.H;
        W_init = zeros(sys.Nt, sys.K);
        for k = 1:sys.K, W_init(:, k) = (channels.h(:, k)' * H_eff_init)'; W_init(:, k) = W_init(:, k)/norm(W_init(:, k)); end
        W_init = W_init * sqrt(sys.Pmax / sys.K);
        
        % --- 1. PROPOSED ---
        W_opt = W_init; Phi_opt = Phi_init;
        for iter = 1:20
            [W_new, ~, err] = Optimize_Beamforming_SCA_Secrecy(sys, channels, W_opt, Phi_opt);
            if ~err, W_opt = W_new; end
            Phi_opt = Optimize_SIM_GA_Secrecy(sys, channels, W_opt, Phi_opt);
        end
        WSR_vs_Gamma_opt_accum(idx) = WSR_vs_Gamma_opt_accum(idx) + Param_cal_Secrecy(sys, channels, W_opt, Phi_opt);
        
        % --- 2. BASELINE (Giữ nguyên Phi_init) ---
        W_rand = W_init;
        [W_new_rand, ~, err_rand] = Optimize_Beamforming_SCA_Secrecy(sys, channels, W_rand, Phi_init);
        if ~err_rand, W_rand = W_new_rand; end
        WSR_vs_Gamma_rand_accum(idx) = WSR_vs_Gamma_rand_accum(idx) + Param_cal_Secrecy(sys, channels, W_rand, Phi_init);
    end
end

WSR_vs_Gamma_opt = WSR_vs_Gamma_opt_accum / sys_base.n_monte;
WSR_vs_Gamma_rand = WSR_vs_Gamma_rand_accum / sys_base.n_monte;

%% 3. LƯU KẾT QUẢ
output_dir = fullfile(pwd, 'output');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end
save_path = fullfile(output_dir, 'Results_P4_Impact_Gamma.mat');
save(save_path, 'Gamma_dBm_list', 'WSR_vs_Gamma_opt', 'WSR_vs_Gamma_rand', 'sys_base');
fprintf('\n--- Đã lưu kết quả Part 4 thành công! ---\n');