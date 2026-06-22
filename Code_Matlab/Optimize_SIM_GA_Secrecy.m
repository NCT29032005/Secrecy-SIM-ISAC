function Phi_new = Optimize_SIM_GA_Secrecy(sys, channels, W_ref, Phi_ref, show_log)
    % Optimized SIM phase shifts using Vectorized Gradient Ascent (A & B Combined)
    if nargin < 5, show_log = false; end
    K = sys.K; L = sys.L; M = sys.M;
    sigma2_u = sys.sigma2_u; sigma2_e = sys.sigma2_e;
    Gamma = sys.Gamma;
    beta_penalty = 1000;
    
    ga_iter = sys.ga_iter; % Lấy số bước nhảy từ config
    Phi_new = Phi_ref;

    for l = 1:L
        % Thực hiện nhiều bước nhảy GA liên tiếp cho lớp l (Phương án A)
        for s = 1:ga_iter
            % 1. Cập nhật ma trận Forward (K) và Backward (L) cho cấu hình hiện tại
            K_mats = cell(L, 1);
            L_mats = cell(L, 1);
            K_mats{1} = channels.H;
            for i = 2:L
                K_mats{i} = channels.Psi * (Phi_new{i-1} * K_mats{i-1});
            end
            L_mats{L} = eye(M);
            for i = L-1:-1:1
                L_mats{i} = L_mats{i+1} * (Phi_new{i+1} * channels.Psi);
            end

            % 2. Tính toán các chỉ số hiện tại (Metrics)
            Theta_curr = eye(M);
            for i = 1:L
                Theta_curr = Phi_new{i} * Theta_curr;
                if i < L, Theta_curr = channels.Psi * Theta_curr; end
            end
            H_eff_curr = Theta_curr * channels.H;
            
            % Tính Secrecy Rate cho tất cả người dùng
            delta_u = zeros(K, K); delta_e = zeros(K, K);
            H_eff_W = H_eff_curr * W_ref; % Precompute
            for k = 1:K
                for j = 1:K
                    delta_u(k, j) = abs(channels.h(:, k)' * H_eff_W(:, j))^2;
                    delta_e(k, j) = abs(channels.v(:, k)' * H_eff_W(:, j))^2;
                end
            end
            I_u = sum(delta_u, 2) - diag(delta_u) + sigma2_u;
            I_e = sum(delta_e, 2) - diag(delta_e) + sigma2_e;
            gamma_u = diag(delta_u) ./ I_u;
            gamma_e = diag(delta_e) ./ I_e;
            R_sec = max(0, log2(1 + gamma_u) - log2(1 + gamma_e));
            [~, k_star] = min(R_sec);
            
            % Tính Sensing Power
            a_H_eff_W = channels.a' * H_eff_W;
            P_sensing = sum(abs(a_H_eff_W).^2);

            % 3. VECTORIZED GRADIENT CALCULATION (Phương án B)
            phi_l = diag(Phi_new{l});
            omega_l = angle(phi_l);
            exp_j_phi = exp(1i * omega_l); % M x 1

            % Precompute common terms for all meta-atoms m
            % [h_k^H Ll]_m -> Vector 1 x M
            h_L_u_kstar = channels.h(:, k_star)' * L_mats{l}; % 1 x M
            h_L_e_kstar = channels.v(:, k_star)' * L_mats{l}; % 1 x M
            a_L = channels.a' * L_mats{l};                   % 1 x M
            
            % [Kl w_j]_m -> Vector M x 1
            K_w = K_mats{l} * W_ref; % M x K
            
            % Vectorized calculation of d_delta for all m
            d_delta_u = zeros(M, K); % M x K (for k_star, each beam j)
            d_delta_e = zeros(M, K);
            
            for j = 1:K
                % UE k_star signal/interference from beam j
                H_eff_W_u = channels.h(:, k_star)' * H_eff_W(:, j);
                d_delta_u(:, j) = 2 * real(1i * (h_L_u_kstar.' .* exp_j_phi .* K_w(:, j)) * conj(H_eff_W_u));
                
                % Eve k_star signal/leakage from beam j
                H_eff_W_e = channels.v(:, k_star)' * H_eff_W(:, j);
                d_delta_e(:, j) = 2 * real(1i * (h_L_e_kstar.' .* exp_j_phi .* K_w(:, j)) * conj(H_eff_W_e));
            end
            
            % Gradient of UE Rate (k_star) - Eq (42) vectorized
            d_I_u = sum(d_delta_u, 2) - d_delta_u(:, k_star);
            d_gamma_u = (d_delta_u(:, k_star) * I_u(k_star) - diag(delta_u(k_star, k_star)) * d_I_u) / (I_u(k_star)^2);
            grad_R_u = (1/(log(2)*(1+gamma_u(k_star)))) * d_gamma_u;
            
            % Gradient of Eve Rate (k_star)
            d_I_e = sum(d_delta_e, 2) - d_delta_e(:, k_star);
            d_gamma_e = (d_delta_e(:, k_star) * I_e(k_star) - diag(delta_e(k_star, k_star)) * d_I_e) / (I_e(k_star)^2);
            grad_R_e = (1/(log(2)*(1+gamma_e(k_star)))) * d_gamma_e;
            
            grad_R_sec_total = grad_R_u - grad_R_e;

            % Vectorized Gradient of Sensing - Eq (45)
            grad_sensing = zeros(M, 1);
            if P_sensing < Gamma
                for k = 1:K
                    a_H_eff_W_k = channels.a' * H_eff_W(:, k);
                    grad_sensing = grad_sensing + 2 * real(1i * (a_L.' .* exp_j_phi .* K_w(:, k)) * conj(a_H_eff_W_k));
                end
            end
            
            % Total Gradient
            grad_total = grad_R_sec_total + beta_penalty * grad_sensing;

            % 4. Backtracking Line Search
            mu = sys.mu_init;
            max_grad = max(abs(grad_total));
            if max_grad > 0, step_dir = grad_total / max_grad; else, step_dir = grad_total; end
            
            Obj_old = R_sec(k_star) + beta_penalty * min(P_sensing - Gamma, 0);
            success = false;
            
            for iter_step = 1:10
                omega_new = omega_l + mu * step_dir;
                Phi_temp_l = diag(exp(1i * omega_new));
                
                % Quick test: Only update Theta for layer l to save time
                Theta_new = eye(M);
                for i = 1:L
                    if i == l, Theta_new = Phi_temp_l * Theta_new;
                    else, Theta_new = Phi_new{i} * Theta_new; end
                    if i < L, Theta_new = channels.Psi * Theta_new; end
                end
                
                [R_sec_new, ~, ~, ~, ~] = quick_metrics(sys, channels, W_ref, Theta_new);
                P_sensing_new = sum(abs(channels.a' * Theta_new * channels.H * W_ref).^2);
                Obj_new = min(R_sec_new) + beta_penalty * min(P_sensing_new - Gamma, 0);
                
                if Obj_new > Obj_old + 1e-6
                    Phi_new{l} = Phi_temp_l;
                    success = true;
                    break;
                end
                mu = mu * 0.5;
            end
            
            if ~success, break; end % Nếu không tìm được bước nhảy tăng, dừng vòng lặp s
        end
        
        if show_log
            % Final check for logging
            Theta_log = eye(M);
            for i = 1:L, Theta_log = Phi_new{i} * Theta_log; if i < L, Theta_log = channels.Psi * Theta_log; end, end
            [R_log, ~, ~, ~, ~] = quick_metrics(sys, channels, W_ref, Theta_log);
            P_log = sum(abs(channels.a' * Theta_log * channels.H * W_ref).^2);
            fprintf('> Lớp %d/%d (ga_iter:%d) | WSR: %.4f | Sens: %.2f\n', l, L, s, min(R_log), P_log);
        end
    end
end

function [R_sec, R_u, R_e, gamma_u, gamma_e] = quick_metrics(sys, channels, W, Theta)
    K = sys.K; sigma2_u = sys.sigma2_u; sigma2_e = sys.sigma2_e;
    H_eff = Theta * channels.H;
    H_eff_W = H_eff * W;
    delta_u = zeros(K, K); delta_e = zeros(K, K);
    for k = 1:K
        for j = 1:K
            delta_u(k, j) = abs(channels.h(:, k)' * H_eff_W(:, j))^2;
            delta_e(k, j) = abs(channels.v(:, k)' * H_eff_W(:, j))^2;
        end
    end
    I_u = sum(delta_u, 2) - diag(delta_u) + sigma2_u;
    I_e = sum(delta_e, 2) - diag(delta_e) + sigma2_e;
    gamma_u = diag(delta_u) ./ I_u;
    gamma_e = diag(delta_e) ./ I_e;
    R_u = log2(1 + gamma_u);
    R_e = log2(1 + gamma_e);
    R_sec = max(0, R_u - R_e);
end
