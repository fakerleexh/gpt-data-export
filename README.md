function analyze_irregular_response()

clc; close all;

    %==================== 配置 ====================%
    filename    = 'irregular.xlsx';
    xlim_hz     = [0, 0.5];      % 频率显示范围（Hz）
    seg_sec     = 60;            % Welch方法每段时长（秒）
    %=============================================%
    
    % 检查文件
    assert(exist(filename, 'file')==2, '未找到文件：%s', filename);
    
    % 读取数据
    T = readtable(filename);
    time        = T.Time;
    roll_noTMD  = T.PtfmRoll;
    pitch_noTMD = T.PtfmPitch;
    roll_TMD    = T.PtfmRoll_TMD;
    pitch_TMD   = T.PtfmPitch_TMD;
    
    % 采样频率
    fs = 1/median(diff(time));
    fprintf('========== 数据信息 ==========\n');
    fprintf('采样频率: %.2f Hz\n', fs);
    fprintf('数据时长: %.2f s\n', time(end)-time(1));
    fprintf('数据点数: %d\n', length(time));
    fprintf('==============================\n\n');
    
    % 去趋势（保留原始数据用于时域图）
    roll_noTMD_raw  = roll_noTMD;
    pitch_noTMD_raw = pitch_noTMD;
    roll_TMD_raw    = roll_TMD;
    pitch_TMD_raw   = pitch_TMD;
    
    roll_noTMD  = detrend(roll_noTMD);
    pitch_noTMD = detrend(pitch_noTMD);
    roll_TMD    = detrend(roll_TMD);
    pitch_TMD   = detrend(pitch_TMD);
    
    %========== 时域对比图 ==========%
    plot_time_domain(time, pitch_noTMD_raw, pitch_TMD_raw, roll_noTMD_raw, roll_TMD_raw);
    
    %========== 分析Pitch ==========%
    analyze_and_plot(pitch_noTMD, pitch_TMD, fs, seg_sec, xlim_hz, 'Pitch');

    fprintf('\n图像已保存：Time_Domain_Comparison.png, Pitch_Analysis.png\n'); % 仅保留Pitch相关提示
end

%================= 时域对比绘图函数 =================%
function plot_time_domain(time, pitch_noTMD, pitch_TMD, roll_noTMD, roll_TMD)
    % 绘制时域对比图
    
    % 计算统计量
    std_pitch_noTMD = std(pitch_noTMD);
    std_pitch_TMD = std(pitch_TMD);
    std_roll_noTMD = std(roll_noTMD);
    std_roll_TMD = std(roll_TMD);
    
    max_pitch_noTMD = max(abs(pitch_noTMD));
    max_pitch_TMD = max(abs(pitch_TMD));
    max_roll_noTMD = max(abs(roll_noTMD));
    max_roll_TMD = max(abs(roll_TMD));
    
    % 计算衰减率
    std_att_pitch = (1 - std_pitch_TMD/std_pitch_noTMD) * 100;
    max_att_pitch = (1 - max_pitch_TMD/max_pitch_noTMD) * 100;
    std_att_roll = (1 - std_roll_TMD/std_roll_noTMD) * 100;
    max_att_roll = (1 - max_roll_TMD/max_roll_noTMD) * 100;
    
    fprintf('========== 时域统计结果 ==========\n');
    fprintf('Pitch:\n');
    fprintf('  标准差 - 无TMD: %.4f deg | 有TMD: %.4f deg | 衰减: %.1f%%\n', ...
        std_pitch_noTMD, std_pitch_TMD, std_att_pitch);
    fprintf('  最大值 - 无TMD: %.4f deg | 有TMD: %.4f deg | 衰减: %.1f%%\n\n', ...
        max_pitch_noTMD, max_pitch_TMD, max_att_pitch);

    fprintf('==================================\n\n');
    
    % 绘图
    figure('Position', [100, 100, 1200, 700], 'Color', 'w');
    
    % 子图1: Pitch时域
    subplot(1,1,1);
    plot(time, pitch_noTMD, 'b-', 'LineWidth', 1.2, 'DisplayName', 'No TMD');
    hold on;
    plot(time, pitch_TMD, 'r-', 'LineWidth', 1.2, 'DisplayName', 'With TMD');
    grid on;
    xlabel('Time (s)', 'FontSize', 11);
    ylabel('Pitch (deg)', 'FontSize', 11);
    title('Pitch 时域响应对比', 'FontSize', 12, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 10);
    
    % 添加文本注释
    ylim_val = ylim;
    text(0.02, 0.95, sprintf('Std减少: %.1f%%\nMax减少: %.1f%%', std_att_pitch, max_att_pitch), ...
        'Units', 'normalized', 'FontSize', 9, 'BackgroundColor', 'white', ...
        'EdgeColor', 'k', 'VerticalAlignment', 'top');
    
    
    % 保存
    saveas(gcf, 'Time_Domain_Comparison.png');
end

%================= 核心分析函数 =================%
function analyze_and_plot(x_noTMD, x_TMD, fs, seg_sec, xlim_hz, dof_name)
    % 对比有无TMD的三种频域分析方法
    
    x_noTMD = x_noTMD(:);
    x_TMD = x_TMD(:);
    N = length(x_noTMD);
    Nfft = 2^nextpow2(N);
    Nseg = max(256, round(fs * seg_sec));
    
    %---------- 方法1: 直接FFT幅值谱 ----------%
    % No TMD
    X = fft(x_noTMD, Nfft);
    X_single = X(1:Nfft/2+1);
    f_fft = (0:Nfft/2)' * fs / Nfft;
    Amp_fft_noTMD = abs(X_single) / N;
    Amp_fft_noTMD(2:end-1) = 2 * Amp_fft_noTMD(2:end-1);
    
    % With TMD
    X = fft(x_TMD, Nfft);
    X_single = X(1:Nfft/2+1);
    Amp_fft_TMD = abs(X_single) / N;
    Amp_fft_TMD(2:end-1) = 2 * Amp_fft_TMD(2:end-1);
    
    %---------- 方法2: PSD (Welch方法) ----------%
    [Pxx_noTMD, f_psd] = pwelch(x_noTMD, hann(Nseg), round(0.5*Nseg), Nfft, fs);
    [Pxx_TMD, ~] = pwelch(x_TMD, hann(Nseg), round(0.5*Nseg), Nfft, fs);
    
    %---------- 方法3: Welch平均幅值谱 ----------%
    [f_welch, Amp_welch_noTMD] = amp_welch(x_noTMD, fs, Nseg, round(0.5*Nseg), hann(Nseg), Nfft);
    [~, Amp_welch_TMD] = amp_welch(x_TMD, fs, Nseg, round(0.5*Nseg), hann(Nseg), Nfft);
    
    %---------- 找主导频率 ----------%
    [peak_fft_noTMD, idx_fft] = max(Amp_fft_noTMD);
    [peak_fft_TMD, ~] = max(Amp_fft_TMD);
    peak_fft_TMD_at_peak = interp1(f_fft, Amp_fft_TMD, f_fft(idx_fft), 'linear');
    
    [peak_psd_noTMD, idx_psd] = max(Pxx_noTMD);
    [peak_psd_TMD, ~] = max(Pxx_TMD);
    peak_psd_TMD_at_peak = interp1(f_psd, Pxx_TMD, f_psd(idx_psd), 'linear');
    
    [peak_welch_noTMD, idx_welch] = max(Amp_welch_noTMD);
    [peak_welch_TMD, ~] = max(Amp_welch_TMD);
    peak_welch_TMD_at_peak = interp1(f_welch, Amp_welch_TMD, f_welch(idx_welch), 'linear');
    
    % 计算衰减率
    att_fft = (1 - peak_fft_TMD_at_peak/peak_fft_noTMD) * 100;
    att_psd = (1 - peak_psd_TMD_at_peak/peak_psd_noTMD) * 100;
    att_welch = (1 - peak_welch_TMD_at_peak/peak_welch_noTMD) * 100;
    
    fprintf('========== %s 频域分析结果 ==========\n', dof_name);
    fprintf('方法1 - FFT幅值谱:\n');
    fprintf('  主导频率: %.4f Hz\n', f_fft(idx_fft));
    fprintf('  无TMD峰值: %.4f deg | 有TMD: %.4f deg | 衰减: %.1f%%\n\n', ...
        peak_fft_noTMD, peak_fft_TMD_at_peak, att_fft);
    
    fprintf('方法2 - PSD (功率谱密度):\n');
    fprintf('  主导频率: %.4f Hz\n', f_psd(idx_psd));
    fprintf('  无TMD峰值: %.4e deg²/Hz | 有TMD: %.4e deg²/Hz | 衰减: %.1f%%\n\n', ...
        peak_psd_noTMD, peak_psd_TMD_at_peak, att_psd);

    fprintf('=====================================\n\n');
    
    %---------- 绘图 ----------%
    figure('Position', [100, 100, 1200, 800], 'Color', 'w');
    
    % 子图1: FFT幅值谱
    subplot(2,1,1);
    plot(f_fft, Amp_fft_noTMD, 'b-', 'LineWidth', 1.8, 'DisplayName', 'No TMD');
    hold on;
    plot(f_fft, Amp_fft_TMD, 'r-', 'LineWidth', 1.8, 'DisplayName', 'With TMD');
    grid on; xlim(xlim_hz);
    xlabel('Frequency (Hz)', 'FontSize', 11); 
    ylabel('Amplitude (deg)', 'FontSize', 11);
    title(sprintf('%s - FFT幅值谱（单边）', dof_name), 'FontSize', 12, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 10);
    % 标注主频
    plot(f_fft(idx_fft), peak_fft_noTMD, 'bo', 'MarkerSize', 7, 'MarkerFaceColor', 'b');
    text(f_fft(idx_fft), peak_fft_noTMD*1.08, sprintf('%.4f Hz', f_fft(idx_fft)), ...
        'HorizontalAlignment', 'center', 'FontSize', 9, 'Color', 'b');
    
    % 子图2: PSD
    subplot(2,1,2);
    plot(f_psd, Pxx_noTMD, 'b-', 'LineWidth', 1.8, 'DisplayName', 'No TMD');
    hold on;
    plot(f_psd, Pxx_TMD, 'r-', 'LineWidth', 1.8, 'DisplayName', 'With TMD');
    grid on; xlim(xlim_hz);
    xlabel('Frequency (Hz)', 'FontSize', 11); 
    ylabel('PSD (deg²/Hz)', 'FontSize', 11);
    title(sprintf('%s - 功率谱密度 (PSD)', dof_name), 'FontSize', 12, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 10);
    % 标注主频
    plot(f_psd(idx_psd), peak_psd_noTMD, 'bo', 'MarkerSize', 7, 'MarkerFaceColor', 'b');
    text(f_psd(idx_psd), peak_psd_noTMD*1.08, sprintf('%.4f Hz', f_psd(idx_psd)), ...
        'HorizontalAlignment', 'center', 'FontSize', 9, 'Color', 'b');
    
    
    % 保存
    saveas(gcf, [dof_name, '_Analysis.png']);
end

%================= Welch平均幅值谱函数 =================%
function [f, Amean] = amp_welch(x, fs, Nseg, nover, win, Nfft)
    % 使用Welch方法计算单边幅值谱（非PSD）
    
    x = x(:);  
    w = win(:);  
    M = numel(w);
    step = M - nover;  
    if step < 1, step = M; end
    
    L = numel(x);
    k = 0;  
    Aacc = 0;
    
    % 幅值标定因子
    scale = sum(w) / 2;
    
    % 分段处理
    for i = 1:step:(L-M+1)
        seg = x(i:i+M-1);
        seg = seg - mean(seg);       % 段内去直流
        X = fft(seg .* w, Nfft);
        X = X(1:Nfft/2+1);
        A = abs(X) / scale;          % 线性幅值谱
        A(2:end-1) = 2 * A(2:end-1); % 单边
        
        if k == 0
            Aacc = A;
        else
            Aacc = Aacc + A;
        end
        k = k + 1;
    end
    
    % 处理数据太短的情况
    if k == 0
        X = fft(x .* w(1:min(M, numel(x))), Nfft);
        X = X(1:Nfft/2+1);
        A = abs(X) / scale;
        A(2:end-1) = 2 * A(2:end-1);
        Aacc = A; 
        k = 1;
    end
    
    Amean = Aacc / k;
    f = (0:Nfft/2)' * fs / Nfft;
end
