function [W_new, t_val, error] = Optimize_Beamforming_SCA_Secrecy(sys, channels, W_ref, Phi_ref)
    % Optimize BS transmit beamforming using SCA for Secrecy SIM-ISAC
    K = sys.K; Nt = sys.Nt; L = sys.L; M = sys.M;
    Pmax = sys.Pmax; Rmin = sys.Rmin; Gamma = sys.Gamma;
    sigma2_u = sys.sigma2_u; sigma2_e = sys.sigma2_e;

    error = false;
    
    %% 1. Compute Effective Channels
    % Theta = Phi_L * Psi * Phi_{L-1} * ... * Psi * Phi_1
    Theta = eye(M);
    for l = 1:L
        Theta = Phi_ref{l} * Theta;
        if l < L
            Theta = channels.Psi * Theta;
        end
    end
    
    H_eff = Theta * channels.H; % M x Nt
    
    h_eff = zeros(Nt, K);
    v_eff = zeros(Nt, K);
    for k = 1:K
        h_eff(:, k) = (channels.h(:, k)' * H_eff)'; % Nt x 1
        v_eff(:, k) = (channels.v(:, k)' * H_eff)'; % Nt x 1
    end
    
    a_eff = (channels.a' * H_eff)'; % Nt x 1 (Sensing effective channel)

    %% 2. SCA Pre-calculations (Reference Point)
    S_u_ref = zeros(K, 1);
    I_u_ref = zeros(K, 1);
    S_e_ref = zeros(K, 1);
    I_e_ref = zeros(K, 1);
    
    for k = 1:K
        S_u_ref(k) = abs(h_eff(:, k)' * W_ref(:, k))^2;
        I_u_ref(k) = sigma2_u;
        for j = 1:K
            if j ~= k
                I_u_ref(k) = I_u_ref(k) + abs(h_eff(:, k)' * W_ref(:, j))^2;
            end
        end
        
        S_e_ref(k) = abs(v_eff(:, k)' * W_ref(:, k))^2;
        I_e_ref(k) = sigma2_e;
        for j = 1:K
            if j ~= k
                I_e_ref(k) = I_e_ref(k) + abs(v_eff(:, k)' * W_ref(:, j))^2;
            end
        end
    end

    %% 3. CVX Optimization
    try
        cvx_begin quiet
            cvx_solver mosek
            variable W(Nt, K) complex
            variable t_sec % Worst-case secrecy rate
            
            expression S_u(K)
            expression I_u(K)
            expression S_e(K)
            expression I_e(K)
            expression R_u_low(K)
            expression R_e_up(K)
            
            % Beampattern gain linearization (Sensing)
            expression P_sensing_low
            P_sensing_ref = 0;
            for k = 1:K
                P_sensing_ref = P_sensing_ref + abs(a_eff' * W_ref(:, k))^2;
            end
            
            P_sensing_low = 0;
            for k = 1:K
                P_sensing_low = P_sensing_low + 2 * real((a_eff' * W_ref(:, k))' * (a_eff' * W(:, k))) - abs(a_eff' * W_ref(:, k))^2;
            end

            for k = 1:K
                
                S_u_linear = 2 * real((h_eff(:, k)' * W_ref(:, k))' * (h_eff(:, k)' * W(:, k))) - abs(h_eff(:, k)' * W_ref(:, k))^2;
                
                % Interference part (convex)
                I_u_val = sigma2_u;
                for j = 1:K
                    if j ~= k
                        I_u_val = I_u_val + square_abs(h_eff(:, k)' * W(:, j));
                    end
                end
                
                % Constant coefficients for log2(1 + x)
                gamma_u_ref = S_u_ref(k) / I_u_ref(k);
                a_k = log2(1 + gamma_u_ref);
                b_k = gamma_u_ref / (log(2) * (1 + gamma_u_ref));
                
                R_u_low(k) = a_k + b_k * ( (S_u_linear / S_u_ref(k)) - (I_u_val / I_u_ref(k))); 
            
                Total_e = sigma2_e;
                for j = 1:K
                    Total_e = Total_e + square_abs(v_eff(:, k)' * W(:, j));
                end
                
               
                log_total_e_up = log2(S_e_ref(k) + I_e_ref(k)) + (1/(log(2)*(S_e_ref(k) + I_e_ref(k)))) * (Total_e - (S_e_ref(k) + I_e_ref(k)));
                
                I_e_linear = sigma2_e;
                for j = 1:K
                    if j ~= k
                        I_e_linear = I_e_linear + 2 * real((v_eff(:, k)' * W_ref(:, j))' * (v_eff(:, k)' * W(:, j))) - abs(v_eff(:, k)' * W_ref(:, j))^2;
                    end
                end
                
                R_e_up(k) = log_total_e_up - (log(I_e_linear) / log(2)); 
            end
            
            maximize t_sec
            subject to
                sum(square_abs(W(:))) <= Pmax;
                for k = 1:K
                    R_u_low(k) - R_e_up(k) >= t_sec;
                    R_u_low(k) >= Rmin;
                end
                P_sensing_low >= Gamma;
                t_sec >= 0;
                
        cvx_end
        
        if contains(cvx_status, 'Infeasible') || contains(cvx_status, 'Failed')
            error = true; W_new = W_ref; t_val = 0;
        else
            W_new = W; t_val = t_sec;
        end
    catch
        error = true; W_new = W_ref; t_val = 0;
    end
end
