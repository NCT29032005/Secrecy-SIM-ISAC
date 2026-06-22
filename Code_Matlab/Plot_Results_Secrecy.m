% Plot_Results_Secrecy.m

clear; clc; close all;

%% 1. CẤU HÌNH DỮ LIỆU & MÀU SẮC
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, 'output'); 
fig_dir    = fullfile(output_dir, 'figure');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

mat_files = dir(fullfile(output_dir, '*.mat'));

if isempty(mat_files)
    warning('Không tìm thấy file .mat nào trong: %s', output_dir);
    return;
end

data = struct();
for i = 1:length(mat_files)
    fprintf('>> Đang load dữ liệu từ: %s\n', mat_files(i).name);
    temp_data = load(fullfile(output_dir, mat_files(i).name));
    fnames = fieldnames(temp_data);
    for j = 1:length(fnames)
        data.(fnames{j}) = temp_data.(fnames{j});
    end
end

% --- Color Palette ---
navy_blue   = [0.00, 0.45, 0.74]; 
dark_green  = [0.47, 0.67, 0.19];  
red_prop    = [0.85, 0.33, 0.10]; 
purple_perf = [0.49, 0.18, 0.56]; 
yellow_gold = [0.93, 0.69, 0.13];
gray_zf     = [0.50, 0.50, 0.50];

% --- Tham số hiển thị ---
mk_size = 8;

%% 2. FIGURE 1: CONVERGENCE ANALYSIS
if isfield(data, 'WSR_conv_opt_final') && isfield(data, 'WSR_conv_rand_final')
    figure('Color','w','Units','pixels','Position', [100 100 560 420]);
    hold on; grid on; box on;
    
    % --- THIẾT LẬP TRỤC ---
    set(gca, 'FontName', 'Palatino Linotype', 'FontSize', 13,'TickLabelInterpreter', 'latex', 'LineWidth', 1);
    set(gca, 'GridLineStyle', '-', 'GridAlpha', 0.3);
    
    y_opt = data.WSR_conv_opt_final;
    y_rand = data.WSR_conv_rand_final;
    x_axis = 1:length(y_opt);
    
    % Proposed SIM-GA Design
    plot(x_axis, y_opt, '-o', 'Color', red_prop, 'LineWidth', 2.0, 'MarkerSize', 6, 'MarkerFaceColor', red_prop,'DisplayName', 'Proposed SIM-GA');
    
    % Random SIM Phase Design
    plot(x_axis, y_rand, '--s', 'Color', navy_blue, 'LineWidth', 2.0, 'MarkerSize', 6, 'MarkerFaceColor', navy_blue, 'DisplayName', 'Random SIM Phase');
        
    xlabel('\textbf{AO Iterations}','Interpreter','latex','FontSize',12,'FontName','Times New Roman');
    ylabel('\textbf{Worst-case Secrecy Rate [bps/Hz]}','Interpreter','latex','FontSize',12,'FontName','Times New Roman');
    
    % Legend
    lgd = legend('show', 'Location', 'SouthEast'); 
    set(lgd, 'Interpreter', 'tex', 'FontSize', 12,'FontName', 'Times New Roman', 'FontWeight', 'normal');
             
    % Save
    print(gcf, fullfile(fig_dir, 'Secrecy_Convergence.pdf'), '-dpdf', '-r300', '-painters');
    fprintf('>> Saved Convergence plot to output/figure/\n');
end

%% 3. FIGURE 2: IMPACT OF SIM LAYERS
if isfield(data, 'L_list') && isfield(data, 'WSR_vs_L_opt')
    figure('Color','w','Units','pixels','Position', [150 150 560 420]);
    hold on; grid on; box on;
    
    % --- THIẾT LẬP TRỤC ---
    set(gca, 'FontName', 'Palatino Linotype', 'FontSize', 13,'TickLabelInterpreter', 'latex', 'LineWidth', 1);
    set(gca, 'GridLineStyle', '-', 'GridAlpha', 0.3);
    
    plot(data.L_list, data.WSR_vs_L_opt, '-o', 'Color', red_prop, 'LineWidth', 2.0, 'MarkerSize', mk_size, 'MarkerFaceColor', red_prop, 'DisplayName', 'Proposed SIM-GA');
    plot(data.L_list, data.WSR_vs_L_rand, '--s', 'Color', navy_blue, 'LineWidth', 2.0, 'MarkerSize', mk_size, 'MarkerFaceColor', navy_blue, 'DisplayName', 'Random SIM Phase');
        
    xlabel('\textbf{Number of SIM Layers ($L$)}','Interpreter','latex','FontSize',12,'FontName','Times New Roman');
    ylabel('\textbf{Avg. Worst-case Secrecy Rate [bps/Hz]}','Interpreter','latex','FontSize',12,'FontName','Times New Roman');
    
    xticks(data.L_list);
    xlim([min(data.L_list), max(data.L_list)]);
    
    % Legend
    lgd = legend('show', 'Location', 'NorthWest'); 
    set(lgd, 'Interpreter', 'tex', 'FontSize', 12,'FontName', 'Times New Roman', 'FontWeight', 'normal');
             
    % Save
    print(gcf, fullfile(fig_dir, 'Secrecy_vs_Layers.pdf'), '-dpdf', '-r300', '-painters');
    fprintf('>> Saved Layers plot to output/figure/\n');
end

%% 4. FIGURE 3: IMPACT OF MAX TRANSMIT POWER (Pmax)
if isfield(data, 'Pmax_dBm_list') && isfield(data, 'WSR_vs_Pmax_opt')
    figure('Color','w','Units','pixels','Position', [200 200 560 420]);
    hold on; grid on; box on;
    
    % --- THIẾT LẬP TRỤC ---
    set(gca, 'FontName', 'Palatino Linotype', 'FontSize', 13,'TickLabelInterpreter', 'latex', 'LineWidth', 1);
    set(gca, 'GridLineStyle', '-', 'GridAlpha', 0.3);
    
    plot(data.Pmax_dBm_list, data.WSR_vs_Pmax_opt, '-o', 'Color', red_prop, 'LineWidth', 2.0,'MarkerSize', mk_size, 'MarkerFaceColor', red_prop,'DisplayName', 'Proposed SIM-GA');
    plot(data.Pmax_dBm_list, data.WSR_vs_Pmax_rand, '--s', 'Color', navy_blue, 'LineWidth', 2.0,'MarkerSize', mk_size, 'MarkerFaceColor', navy_blue, 'DisplayName', 'Random SIM Phase');
        
    xlabel('\textbf{Max Transmit Power $P_{max}$ [dBm]}','Interpreter','latex','FontSize',12,'FontName','Times New Roman');
    ylabel('\textbf{Avg. Worst-case Secrecy Rate [bps/Hz]}','Interpreter','latex','FontSize',12,'FontName','Times New Roman');
    
    xticks(data.Pmax_dBm_list);
    
    % Legend
    lgd = legend('show', 'Location', 'NorthWest'); 
    set(lgd, 'Interpreter', 'tex', 'FontSize', 12,'FontName', 'Times New Roman', 'FontWeight', 'normal');
             
    % Save
    print(gcf, fullfile(fig_dir, 'Secrecy_vs_Pmax.pdf'), '-dpdf', '-r300', '-painters');
    fprintf('>> Saved Pmax plot to output/figure/\n');
end

%% 5. FIGURE 4: IMPACT OF SENSING THRESHOLD (Gamma)
if isfield(data, 'Gamma_dBm_list') && isfield(data, 'WSR_vs_Gamma_opt') && isfield(data, 'WSR_vs_Gamma_rand')
    figure('Color','w','Units','pixels','Position', [250 250 560 420]);
    hold on; grid on; box on;
    
    set(gca, 'FontName', 'Palatino Linotype', 'FontSize', 13,'TickLabelInterpreter', 'latex', 'LineWidth', 1);
    set(gca, 'GridLineStyle', '-', 'GridAlpha', 0.3);
    
    % --- VẼ ĐỒ THỊ ---
    % Proposed
    plot(data.Gamma_dBm_list, data.WSR_vs_Gamma_opt, '-o', 'Color', red_prop, 'LineWidth', 2.0,'MarkerSize', mk_size, 'MarkerFaceColor', red_prop, 'DisplayName', 'Proposed SIM-GA');
    
    % Baseline
    plot(data.Gamma_dBm_list, data.WSR_vs_Gamma_rand, '--s', 'Color', navy_blue, 'LineWidth', 2.0,'MarkerSize', mk_size, 'MarkerFaceColor', navy_blue, 'DisplayName', 'Random SIM Phase');
        
    % --- NHÃN VÀ LEGEND ---
    xlabel('\textbf{Sensing Beampattern Gain Threshold $\Gamma$ [dBm]}','Interpreter','latex','FontSize',12,'FontName','Times New Roman');
    ylabel('\textbf{Avg. Worst-case Secrecy Rate [bps/Hz]}','Interpreter','latex','FontSize',12,'FontName','Times New Roman');
    
    xticks(data.Gamma_dBm_list);
    legend('Location', 'southwest', 'Interpreter', 'latex', 'FontSize', 11);
   
    % Save
    fig_path = fullfile(fig_dir, 'Secrecy_vs_Gamma.pdf');
    print(gcf, fig_path, '-dpdf', '-r300', '-painters');
    fprintf('>> Saved Gamma plot to %s\n', fig_path);
end
fprintf('>> All plots updated and saved in %s\n', fig_dir);