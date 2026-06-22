% =========================================================================
% File: Main_P3_Impact_Pmax.m
% Mục đích: Khảo sát ảnh hưởng của Công suất phát (Pmax: 20 -> 40 dBm)
%           So sánh Proposed (tối ưu W, Phi) và Baseline (tối ưu W, Phi ngẫu nhiên)
% =========================================================================
clear; clc; close all;
%% 1. KHỞI TẠO CẤU HÌNH GỐC
SIM_ISAC_config_Secrecy; 
sys_base = sys; 
fprintf('\n>>> RUNNING PART 3: IMPACT OF Pmax (%d Monte Carlo) <<<\n', sys.n_monte);

%% 2. VÒNG LẶP MONTE CARLO
Pmax_dBm_list = 20:5:40; 
WSR_vs_Pmax_opt_accum = zeros(length(Pmax_dBm_list), 1);
WSR_vs_Pmax_rand_accum = zeros(length(Pmax_dBm_list), 1); 

for m = 1:sys_base.n_monte
    fprintf('Part 3 - Monte Carlo [%d / %d]\n', m, sys_base.n_monte);
    [channels, ~] = genChannel_Secrecy(sys_base); % Cố định kênh cho lượt chạy này
    
    % Khởi tạo pha SIM ngẫu nhiên dùng chung cho cả 2 đường
    Phi_init = cell(sys_base.L, 1);
    for l = 1:sys_base.L, Phi_init{l} = diag(exp(1i * 2 * pi * rand(sys_base.M, 1))); end
    
    for idx = 1:length(Pmax_dBm_list)
        sys = sys_base; 
        sys.Pmax = 10^(Pmax_dBm_list(idx) / 10);
        
        % Khởi tạo W theo MRT (Maximum Ratio Transmission) cho công suất hiện tại
        Theta_init = eye(sys.M);
        for l = 1:sys.L, Theta_init = Phi_init{l} * Theta_init; if l < sys.L, Theta_init = channels.Psi * Theta_init; end, end
        H_eff_init = Theta_init * channels.H;
        W_init = zeros(sys.Nt, sys.K);
        for k = 1:sys.K, W_init(:, k) = (channels.h(:, k)' * H_eff_init)'; W_init(:, k) = W_init(:, k)/norm(W_init(:, k)); end
        W_init = W_init * sqrt(sys.Pmax / sys.K);
        
        % --- 1. TỐI ƯU PROPOSED ---
        W_opt = W_init; Phi_opt = Phi_init;
        for iter = 1:15 
            [W_new, ~, err] = Optimize_Beamforming_SCA_Secrecy(sys, channels, W_opt, Phi_opt);
            if ~err, W_opt = W_new; end
            Phi_opt = Optimize_SIM_GA_Secrecy(sys, channels, W_opt, Phi_opt);
        end
        WSR_vs_Pmax_opt_accum(idx) = WSR_vs_Pmax_opt_accum(idx) + Param_cal_Secrecy(sys, channels, W_opt, Phi_opt);
        
        % --- 2. BASELINE (Giữ nguyên Phi_init) ---
        W_rand = W_init;
        [W_new_rand, ~, ~] = Optimize_Beamforming_SCA_Secrecy(sys, channels, W_rand, Phi_init);
        if ~isempty(W_new_rand) && ~islogical(W_new_rand), W_rand = W_new_rand; end
        WSR_vs_Pmax_rand_accum(idx) = WSR_vs_Pmax_rand_accum(idx) + Param_cal_Secrecy(sys, channels, W_rand, Phi_init);
    end
end

WSR_vs_Pmax_opt = WSR_vs_Pmax_opt_accum / sys_base.n_monte;
WSR_vs_Pmax_rand = WSR_vs_Pmax_rand_accum / sys_base.n_monte;

%% 3. LƯU KẾT QUẢ
output_dir = fullfile(pwd, 'output');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end
save_path = fullfile(output_dir, 'Results_P3_Impact_Pmax.mat');
save(save_path, 'Pmax_dBm_list', 'WSR_vs_Pmax_opt', 'WSR_vs_Pmax_rand', 'sys_base');
fprintf('\n--- Đã lưu kết quả Part 3 thành công! ---\n');