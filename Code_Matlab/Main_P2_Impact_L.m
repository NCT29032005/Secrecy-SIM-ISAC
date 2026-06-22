% =========================================================================
% File: Main_P2_Impact_L.m
% Mục đích: Khảo sát ảnh hưởng của số lớp SIM (L = 1 đến 6)
%           So sánh giữa Thuật toán Proposed và Random Baseline
% =========================================================================
clear; clc; close all;
%% 1. KHỞI TẠO CẤU HÌNH GỐC
SIM_ISAC_config_Secrecy; 
sys_base = sys; 
fprintf('\n>>> RUNNING PART 2: IMPACT OF SIM LAYERS (%d Monte Carlo) <<<\n', sys.n_monte);

%% 2. VÒNG LẶP MONTE CARLO
L_list = [1, 2, 3, 4, 5, 6];
WSR_vs_L_opt_accum = zeros(length(L_list), 1);
WSR_vs_L_rand_accum = zeros(length(L_list), 1); % Accumulator cho đường Random

for m = 1:sys_base.n_monte
    fprintf('Part 2 - Monte Carlo [%d / %d]\n', m, sys_base.n_monte);
    
    for idx = 1:length(L_list)
        sys = sys_base; L_val = L_list(idx);
        sys.L = L_val; sys.d_layer = sys.T_total / L_val; 
        
        [channels, ~] = genChannel_Secrecy(sys); 
        
        % Khởi tạo pha SIM ngẫu nhiên (dùng chung cho cả 2 trường hợp)
        Phi_init = cell(L_val, 1);
        for l = 1:L_val, Phi_init{l} = diag(exp(1i * 2 * pi * rand(sys.M, 1))); end
        
        % Khởi tạo Beamforming MRT (dùng chung)
        Theta_init = eye(sys.M);
        for l = 1:L_val, Theta_init = Phi_init{l} * Theta_init; if l < L_val, Theta_init = channels.Psi * Theta_init; end, end
        H_eff_init = Theta_init * channels.H;
        W_init = zeros(sys.Nt, sys.K);
        for k = 1:sys.K, W_init(:, k) = (channels.h(:, k)' * H_eff_init)'; W_init(:, k) = W_init(:, k)/norm(W_init(:, k)); end
        W_init = W_init * sqrt(sys.Pmax / sys.K);
        
        % --- 1. TỐI ƯU PROPOSED (W và Phi) ---
        W_opt = W_init; Phi_opt = Phi_init;
        for iter = 1:15
            [W_new_opt, ~, err_W] = Optimize_Beamforming_SCA_Secrecy(sys, channels, W_opt, Phi_opt);
            if ~err_W, W_opt = W_new_opt; end
            Phi_opt = Optimize_SIM_GA_Secrecy(sys, channels, W_opt, Phi_opt);
        end
        WSR_vs_L_opt_accum(idx) = WSR_vs_L_opt_accum(idx) + Param_cal_Secrecy(sys, channels, W_opt, Phi_opt);
        
        % --- 2. BASELINE (Chỉ tối ưu W, giữ nguyên Phi_init) ---
        W_rand = W_init;
        [W_new_rand, ~, err_W_rand] = Optimize_Beamforming_SCA_Secrecy(sys, channels, W_rand, Phi_init);
        if ~err_W_rand, W_rand = W_new_rand; end
        WSR_vs_L_rand_accum(idx) = WSR_vs_L_rand_accum(idx) + Param_cal_Secrecy(sys, channels, W_rand, Phi_init);
    end
end

WSR_vs_L_opt = WSR_vs_L_opt_accum / sys_base.n_monte;
WSR_vs_L_rand = WSR_vs_L_rand_accum / sys_base.n_monte;

%% 3. LƯU KẾT QUẢ
output_dir = fullfile(pwd, 'output');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end
save_path = fullfile(output_dir, 'Results_P2_Impact_L.mat');
save(save_path, 'L_list', 'WSR_vs_L_opt', 'WSR_vs_L_rand', 'sys_base');
fprintf('\n--- Đã lưu kết quả Part 2 thành công! ---\n');