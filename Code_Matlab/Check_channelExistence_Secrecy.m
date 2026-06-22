function [channels] = Check_channelExistence_Secrecy(sys)
    % Check if channel file exists, if not, create and save it.
    % This ensures consistent channels across different simulation runs.
    
    Nt = sys.Nt; M = sys.M; L = sys.L; K = sys.K;
    Monte = sys.n_monte;
    
    data_dir = './output/data_secrecy/';
    file_name = sprintf('%schannel_K%d_Nt%d_M%d_L%d_Monte%d.mat', data_dir, K, Nt, M, L, Monte);
    
    if ~exist(data_dir, 'dir')
        mkdir(data_dir);
    end

    if exist(file_name, 'file') == 0
        fprintf('Không tìm thấy tệp kênh. Đang tạo các kênh mới cho %d lần thử nghiệm Monte Carlo...\n', Monte);
        
        all_channels = cell(Monte, 1);
        for m = 1:Monte
            trial_channels = genChannel_Secrecy(sys);   % Sửa: chỉ 1 output
            all_channels{m} = trial_channels;
        end
        
        save(file_name, 'all_channels', 'sys');   % Lưu thêm sys để kiểm tra sau
        fprintf('Channels saved to: %s\n', file_name);
    else
        fprintf('Đã tìm thấy tệp kênh. Đang tải các kênh hiện có từ: %s\n', file_name);
        data = load(file_name);
        all_channels = data.all_channels;
    end
    
    channels = all_channels;
end
