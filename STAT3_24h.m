
%% TDT 24-hour intermittent fiber photometry dF/F analysis
% Only export summary, no full trace saving

clear; clc; close all;

%% ================= USER SETTINGS =================
SDKPATH   = '/Users/jashinchan/Desktop/Research/TDTMatlabSDK';       % CHANGE THIS
BLOCKPATH = '/Volumes/Koharu/Core_Research/STAT3/260701_4-193_4-194_1e13_LPS_Post';     % CHANGE THIS
OUTDIR    = '/Volumes/Koharu/Core_Research/STAT3/Analysis';    % CHANGE THIS

on_duration_sec  = 5 * 60;      % 5 min recording
cycle_period_sec = 15 * 60;     % 5 min on + 10 min off
n_cycles = 96;

summary_start_sec = 60;         % middle 3 min start
summary_end_sec   = 240;        % middle 3 min end
smooth_filter_order = 1000;     % matches preprocess_sleep_data.m smoothing; no downsampling
zscore_reference_cycles = 1:4;  % fixed first-hour reference; uses middle 3 min of each cycle

channels = {
    'A', 'x405A', 'x465A';
    'C', 'x405C', 'x465C'
};

excel_file = fullfile(OUTDIR, 'TDT_dFF_summary.xlsx');
mat_file   = fullfile(OUTDIR, 'TDT_dFF_summary.mat');

%% ================= SETUP =================
if ~exist(OUTDIR, 'dir')
    mkdir(OUTDIR);
end

addpath(genpath(SDKPATH));

if exist(excel_file, 'file')
    delete(excel_file);
end

fprintf('Importing TDT block...\n');
data = TDTbin2mat(BLOCKPATH);
fprintf('Import complete.\n');

all_summary = table();

%% ================= MAIN ANALYSIS =================
for ch = 1:size(channels,1)

    label    = channels{ch,1};
    store405 = channels{ch,2};
    store465 = channels{ch,3};

    fprintf('\nProcessing channel %s: %s + %s\n', label, store405, store465);

    if ~isfield(data.streams, store405)
        warning('Stream %s not found. Skipping channel %s.', store405, label);
        continue
    end

    if ~isfield(data.streams, store465)
        warning('Stream %s not found. Skipping channel %s.', store465, label);
        continue
    end

    F405 = double(data.streams.(store405).data(:));
    F465 = double(data.streams.(store465).data(:));
    fs   = double(data.streams.(store465).fs);

    n = min(length(F405), length(F465));
    F405 = F405(1:n);
    F465 = F465(1:n);

    t = (0:n-1)' / fs;

    cycle_results = struct( ...
        'Cycle', {}, ...
        'Clock_time_hr', {}, ...
        'Window_start_s', {}, ...
        'Window_end_s', {}, ...
        'Summary_start_s', {}, ...
        'Summary_end_s', {}, ...
        'Fs_Hz', {}, ...
        'Mean_405_middle3min', {}, ...
        'Mean_465_middle3min', {}, ...
        'Mean_fitted405_middle3min', {}, ...
        'DFF_summary', {}, ...
        'Is_zscore_reference_cycle', {}, ...
        'Fit_slope', {}, ...
        'Fit_intercept', {});
    channel_reference_dff_for_zscore = [];
    zscore_reference_label = sprintf('cycles_%d_to_%d_middle3min', ...
        zscore_reference_cycles(1), zscore_reference_cycles(end));

    for cyc = 1:n_cycles

        t_start = (cyc-1) * cycle_period_sec;
        t_end   = t_start + on_duration_sec;

        idx = t >= t_start & t < t_end;

        if sum(idx) < fs * 10
            fprintf('Skipping cycle %d channel %s: not enough data.\n', cyc, label);
            continue
        end

        tc = t(idx) - t_start;

        f405 = F405(idx);
        f465 = F465(idx);

        %% ===== Fit 405 control to 465 signal =====
        valid_fit = isfinite(f405) & isfinite(f465);

        if sum(valid_fit) < 10
            fprintf('Skipping cycle %d channel %s: not enough valid points.\n', cyc, label);
            continue
        end

        p = polyfit(f405(valid_fit), f465(valid_fit), 1);
        fitted405 = polyval(p, f405);

        %% ===== Smooth fitted 405 control =====
        fitted405 = smooth_trace_zero_phase(fitted405, smooth_filter_order);

        %% ===== Calculate dF/F =====
        dff_percent = (f465 - fitted405) ./ fitted405 * 100;

        %% ===== Smooth dF/F =====
        dff_percent = smooth_trace_zero_phase(dff_percent, smooth_filter_order);

        %% ===== Middle 3 min summary only =====
        summary_idx = tc >= summary_start_sec & tc < summary_end_sec;

        dff_summary = dff_percent(summary_idx);
        f405_summary = f405(summary_idx);
        f465_summary = f465(summary_idx);
        fitted405_summary = fitted405(summary_idx);

        valid_summary = isfinite(dff_summary);

        if sum(valid_summary) < 10
            fprintf('Skipping summary cycle %d channel %s: not enough valid dF/F.\n', cyc, label);
            continue
        end

        dff_summary = dff_summary(valid_summary);
        f405_summary = f405_summary(valid_summary);
        f465_summary = f465_summary(valid_summary);
        fitted405_summary = fitted405_summary(valid_summary);

        is_zscore_reference_cycle = ismember(cyc, zscore_reference_cycles);

        if is_zscore_reference_cycle
            channel_reference_dff_for_zscore = [channel_reference_dff_for_zscore; dff_summary];
        end

        mean_405 = mean(f405_summary, 'omitnan');
        mean_465 = mean(f465_summary, 'omitnan');
        mean_fitted405 = mean(fitted405_summary, 'omitnan');

        clock_time_hr = t_start / 3600;

        cycle_results(end+1).Cycle = cyc;
        cycle_results(end).Clock_time_hr = clock_time_hr;
        cycle_results(end).Window_start_s = t_start;
        cycle_results(end).Window_end_s = t_end;
        cycle_results(end).Summary_start_s = summary_start_sec;
        cycle_results(end).Summary_end_s = summary_end_sec;
        cycle_results(end).Fs_Hz = fs;
        cycle_results(end).Mean_405_middle3min = mean_405;
        cycle_results(end).Mean_465_middle3min = mean_465;
        cycle_results(end).Mean_fitted405_middle3min = mean_fitted405;
        cycle_results(end).DFF_summary = dff_summary;
        cycle_results(end).Is_zscore_reference_cycle = is_zscore_reference_cycle;
        cycle_results(end).Fit_slope = p(1);
        cycle_results(end).Fit_intercept = p(2);

    end

    channel_summary = table();

    zscore_mean_dff = mean(channel_reference_dff_for_zscore, 'omitnan');
    zscore_std_dff = std(channel_reference_dff_for_zscore, 'omitnan');

    if isempty(channel_reference_dff_for_zscore) || ~isfinite(zscore_std_dff) || zscore_std_dff == 0
        zscore_mean_dff = NaN;
        zscore_std_dff = NaN;
        warning('Channel %s has no usable first-hour reference dF/F variability for z-scoring.', label);
    end

    for result_idx = 1:numel(cycle_results)

        dff_summary = cycle_results(result_idx).DFF_summary;

        %% ===== Summary parameters =====
        mean_dff   = mean(dff_summary, 'omitnan');
        median_dff = median(dff_summary, 'omitnan');
        max_dff    = max(dff_summary);
        min_dff    = min(dff_summary);
        std_dff    = std(dff_summary, 'omitnan');

        z_dff_summary = (dff_summary - zscore_mean_dff) ./ zscore_std_dff;

        mean_z_dff   = mean(z_dff_summary, 'omitnan');
        median_z_dff = median(z_dff_summary, 'omitnan');
        max_z_dff    = max(z_dff_summary);
        min_z_dff    = min(z_dff_summary);
        std_z_dff    = std(z_dff_summary, 'omitnan');

        temp_summary = table( ...
            string(label), cycle_results(result_idx).Cycle, cycle_results(result_idx).Clock_time_hr, ...
            cycle_results(result_idx).Window_start_s, cycle_results(result_idx).Window_end_s, ...
            cycle_results(result_idx).Summary_start_s, cycle_results(result_idx).Summary_end_s, ...
            cycle_results(result_idx).Fs_Hz, smooth_filter_order, ...
            cycle_results(result_idx).Mean_405_middle3min, ...
            cycle_results(result_idx).Mean_465_middle3min, ...
            cycle_results(result_idx).Mean_fitted405_middle3min, ...
            mean_dff, median_dff, max_dff, min_dff, std_dff, ...
            mean_z_dff, median_z_dff, max_z_dff, min_z_dff, std_z_dff, ...
            zscore_mean_dff, zscore_std_dff, string(zscore_reference_label), ...
            cycle_results(result_idx).Is_zscore_reference_cycle, ...
            cycle_results(result_idx).Fit_slope, cycle_results(result_idx).Fit_intercept, ...
            'VariableNames', { ...
            'Channel','Cycle','Clock_time_hr', ...
            'Window_start_s','Window_end_s', ...
            'Summary_start_s','Summary_end_s', ...
            'Fs_Hz','Smooth_filter_order', ...
            'Mean_405_middle3min','Mean_465_middle3min','Mean_fitted405_middle3min', ...
            'Mean_dFF_percent_middle3min','Median_dFF_percent_middle3min', ...
            'Max_dFF_percent_middle3min','Min_dFF_percent_middle3min', ...
            'STD_dFF_percent_middle3min', ...
            'Mean_z_dFF_middle3min','Median_z_dFF_middle3min', ...
            'Max_z_dFF_middle3min','Min_z_dFF_middle3min', ...
            'STD_z_dFF_middle3min', ...
            'Zscore_reference_mean_dFF_percent','Zscore_reference_STD_dFF_percent', ...
            'Zscore_reference_label','Is_zscore_reference_cycle', ...
            'Fit_slope','Fit_intercept'});

        channel_summary = [channel_summary; temp_summary];

    end

    all_summary = [all_summary; channel_summary];

    writetable(channel_summary, excel_file, 'Sheet', ['Summary_Channel_' label]);

end

%% ================= EXPORT =================
writetable(all_summary, excel_file, 'Sheet', 'Summary_All');

save(mat_file, 'all_summary');

fprintf('\nAnalysis complete.\n');
fprintf('Summary Excel saved to:\n%s\n', excel_file);
fprintf('Summary MAT saved to:\n%s\n', mat_file);

function smoothed_trace = smooth_trace_zero_phase(trace, filter_order)
    smoothed_trace = double(trace(:));

    if filter_order <= 1
        return
    end

    if numel(smoothed_trace) <= 3 * filter_order
        warning('Skipping smoothing because the trace is too short for filter order %d.', filter_order);
        return
    end

    if any(~isfinite(smoothed_trace))
        warning('Skipping smoothing because the trace contains non-finite values.');
        return
    end

    mean_filter = ones(filter_order, 1) / filter_order;
    smoothed_trace = filtfilt(mean_filter, 1, smoothed_trace);
end
