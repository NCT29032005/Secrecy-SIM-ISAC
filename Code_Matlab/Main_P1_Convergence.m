% =========================================================================
% File: Main_P1_Convergence.m
% Mục đích: Khảo sát sự hội tụ của thuật toán (WSR vs. Iterations)
% =========================================================================
clear; clc; close all;

%% 1. KHỞI TẠO CẤU HÌNH GỐC
SIM_ISAC_config_Secrecy; 
fprintf('\n>>> RUNNING PART 1: CONVERGENCE ANALYSIS (%d Monte Carlo) <<<\n', sys.n_monte);

%% 2. TẠO TẬP HỢP KÊNH TRUYỀN 
all_channels_conv = Check_channelExistence_Secrecy(sys); 

%% 3. VÒNG LẶP MONTE CARLO
WSR_conv_opt_accum = zeros(sys.itermax, 1); 
WSR_conv_rand_accum = zeros(sys.itermax, 1);
max_iter_reached = 0;
    
for m = 1:sys.n_monte
    fprintf('Part 1 - Monte Carlo [%d / %d]\n', m, sys.n_monte);
    channels = all_channels_conv{m};
    
    % Khởi tạo pha ngẫu nhiên cho SIM
    Phi_init = cell(sys.L, 1);                                                                                                                                   
    for l = 1:sys.L, Phi_init{l} = diag(exp(1i * 0.05 * randn(sys.M, 1))); end
    
    % Tính kênh hiệu dụng ban đầu
    Theta_init = eye(sys.M);
    for l = 1:sys.L
        Theta_init = Phi_init{l} * Theta_init; 
        if l < sys.L, Theta_init = channels.Psi * Theta_init; end
    end
    H_eff_init = Theta_init * channels.H;
    
    % Khởi tạo MRT cho ma trận W (Cực kỳ quan trọng để tránh Infeasible)
    W_init = zeros(sys.Nt, sys.K);
    for k = 1:sys.K
        h_k_eff = (channels.h(:, k)' * H_eff_init)'; 
        W_init(:, k) = h_k_eff / norm(h_k_eff); 
    end
    W_init = W_init * sqrt(sys.Pmax / sys.K);
    
    W_ref_opt = W_init; Phi_ref_opt = Phi_init;
    W_ref_rand = W_init;
    
    WSR_temp_opt = zeros(sys.itermax, 1);
    WSR_temp_rand = zeros(sys.itermax, 1);
    
    for iter = 1:sys.itermax
        % Thuật toán đề xuất (Tối ưu cả W và Phi)
        [W_new_opt, ~, err_W] = Optimize_Beamforming_SCA_Secrecy(sys, channels, W_ref_opt, Phi_ref_opt);
        if ~err_W, W_ref_opt = W_new_opt; end
        Phi_ref_opt = Optimize_SIM_GA_Secrecy(sys, channels, W_ref_opt, Phi_ref_opt, false); 
        
        % Baseline: Chỉ tối ưu W, giữ nguyên Phi ngẫu nhiên ban đầu
        [W_new_rand, ~, ~] = Optimize_Beamforming_SCA_Secrecy(sys, channels, W_ref_rand, Phi_init);
        if ~isempty(W_new_rand) && ~islogical(W_new_rand), W_ref_rand = W_new_rand; end
        
        WSR_temp_opt(iter) = Param_cal_Secrecy(sys, channels, W_ref_opt, Phi_ref_opt);
        WSR_temp_rand(iter) = Param_cal_Secrecy(sys, channels, W_ref_rand, Phi_init);
        
        % Điều kiện dừng sớm
        if iter > 1 && abs(WSR_temp_opt(iter) - WSR_temp_opt(iter-1)) < sys.tol
            WSR_temp_opt(iter:end) = WSR_temp_opt(iter);
            WSR_temp_rand(iter:end) = WSR_temp_rand(iter);
            max_iter_reached = max(max_iter_reached, iter);
            break; 
        end
        max_iter_reached = max(max_iter_reached, iter);
    end
    WSR_conv_opt_accum = WSR_conv_opt_accum + WSR_temp_opt;
    WSR_conv_rand_accum = WSR_conv_rand_accum + WSR_temp_rand;
end

WSR_conv_opt_final = WSR_conv_opt_accum(1:max_iter_reached) / sys.n_monte;
WSR_conv_rand_final = WSR_conv_rand_accum(1:max_iter_reached) / sys.n_monte;

%% 4. LƯU KẾT QUẢ ĐỘC LẬP
output_dir = fullfile(pwd, 'output');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end
save_path = fullfile(output_dir, 'Results_P1_Convergence.mat');
save(save_path, 'WSR_conv_opt_final', 'WSR_conv_rand_final', 'sys');
fprintf('\n--- Đã lưu kết quả Part 1 thành công tại: %s ---\n', save_path);