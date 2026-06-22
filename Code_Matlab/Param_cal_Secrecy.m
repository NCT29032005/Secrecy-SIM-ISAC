function [WSR, R_sec, P_trans, P_pa, P_SIM] = Param_cal_Secrecy(sys, channels, W, Phi)
% Param_cal_Secrecy: Calculates Secrecy metrics for ISAC SIM system
% Based on the structure of Simulations/Param_cal.m

K = sys.K; 
L = sys.L;
M = sys.M;
sigma2_u = sys.sigma2_u; 
sigma2_e = sys.sigma2_e;

% Calculate Effective Channel Matrix G (Theta)
Theta = eye(M);
for l = 1:L
    Theta = Phi{l} * Theta;
    if l < L
        Theta = channels.Psi * Theta;
    end
end
H_eff = Theta * channels.H;

% --- Secrecy Rate Calculation ---
R_sec = zeros(K, 1);
for k = 1:K
    % Signal and Interference at User k
    h_k_eff = channels.h(:, k)' * H_eff;
    S_u = abs(h_k_eff * W(:, k))^2;
    I_u = sigma2_u;
    for j = 1:K
        if j ~= k
            I_u = I_u + abs(h_k_eff * W(:, j))^2;
        end
    end
    
    % Signal and Interference at Eavesdropper k
    v_k_eff = channels.v(:, k)' * H_eff;
    S_e = abs(v_k_eff * W(:, k))^2;
    I_e = sigma2_e;
    for j = 1:K
        if j ~= k
            I_e = I_e + abs(v_k_eff * W(:, j))^2;
        end
    end
    
    % Secrecy Rate for User k
    R_sec(k) = max(0, log2(1 + S_u/I_u) - log2(1 + S_e/I_e));
end

WSR = min(R_sec); % Worst-case Secrecy Rate

% --- Power Calculation (Optional, based on Param_cal.m structure) ---
% Note: Using parameters from sys if they exist, otherwise providing defaults or placeholders
if isfield(sys, 'rho'), rho = sys.rho; else, rho = 0.35; end
if isfield(sys, 'P_meta'), P_meta = sys.P_meta; else, P_meta = 0; end
if isfield(sys, 'P_ctrl'), P_ctrl = sys.P_ctrl; else, P_ctrl = 0; end
if isfield(sys, 'P0'), P0 = sys.P0; else, P0 = 0; end

P_pa = (1 / rho) * norm(W, 'fro')^2;
P_SIM = L * M * P_meta + P_ctrl;
P_trans = P_pa + P_SIM + P0;

end
