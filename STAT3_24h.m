%% TDT 24-hour intermittent fiber photometry dF/F analysis
% Export cycle-level summary tables; data files and full traces are not saved.

clear; clc; close all;

%% ================= USER SETTINGS =================
SDKPATH   = '/Users/jashinchan/Desktop/Research/TDTMatlabSDK';       % CHANGE THIS
BLOCKPATH = '/Volumes/Koharu/Core_Research/STAT3/260701_4-193_4-194_1e13_LPS_Post';     % CHANGE THIS
OUTDIR    = '/Volumes/Koharu/Core_Research/STAT3/Analysis';          % CHANGE THIS

on_duration_sec  = 5 * 60;      % 5 min recording
cycle_period_sec = 15 * 60;     % 5 min on + 10 min off
n_cycles = 96;

summary_start_sec = 60;         % middle 3 min start
summary_end_sec   = 240;        % middle 3 min end
smooth_filter_sec = 1;          % reference-style moving average, defined in seconds

zscore_reference_cycles = 1:4;  % fixed first-hour reference; uses middle 3 min of each cycle

% Use a fixed first-hour fit by default so sustained 24 h signal changes are
% not absorbed by a new per-cycle intercept. Set to 'per_cycle' for the old
% behavior.
fit_mode = 'reference_cycles';  % 'reference_cycles' or 'per_cycle'
fit_reference_cycles = zscore_reference_cycles;
fit_reference_start_sec = summary_start_sec;
fit_reference_end_sec = summary_end_sec;

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
    fs405 = double(data.streams.(store405).fs);
    fs465 = double(data.streams.(store465).fs);

    if abs(fs405 - fs465) > max(eps(fs465), fs465 * 1e-6)
        warning('Channel %s has mismatched fs values: %s=%g Hz, %s=%g Hz. Skipping.', ...
            label, store405, fs405, store465, fs465);
        continue
    end

    fs = fs465;

    if abs(length(F405) - length(F465)) > fs
        warning('Channel %s stream lengths differ by more than 1 s: %s=%d, %s=%d. Truncating to shorter stream.', ...
            label, store405, length(F405), store465, length(F465));
    end

    n = min(length(F405), length(F465));
    F405 = F405(1:n);
    F465 = F465(1:n);

    t = (0:n-1)' / fs;

    expected_total_sec = n_cycles * cycle_period_sec;
    expected_last_window_end_sec = (n_cycles - 1) * cycle_period_sec + on_duration_sec;
    stream_duration_sec = (n - 1) / fs;

    if stream_duration_sec < expected_last_window_end_sec
        warning(['Channel %s stream duration is %.2f h, shorter than the last expected active window end %.2f h. ', ...
            'Later cycles will be missing or misaligned if the block is not continuous.'], ...
            label, stream_duration_sec / 3600, expected_last_window_end_sec / 3600);
    elseif stream_duration_sec < expected_total_sec * 0.95
        warning(['Channel %s stream duration is %.2f h, shorter than the expected 24 h cycle grid %.2f h. ', ...
            'Verify whether the TDT block is continuous or on-period snippets only.'], ...
            label, stream_duration_sec / 3600, expected_total_sec / 3600);
    end

    smooth_filter_order = max(1, round(smooth_filter_sec * fs));

    [reference_fit, reference_fit_ok] = compute_reference_fit( ...
        F405, F465, t, fs, fit_reference_cycles, cycle_period_sec, ...
        fit_reference_start_sec, fit_reference_end_sec);

    fit_mode_used = fit_mode;
    if strcmpi(fit_mode, 'reference_cycles') && ~reference_fit_ok
        warning('Channel %s could not build a reference-cycle fit. Falling back to per-cycle fits.', label);
        fit_mode_used = 'per_cycle_fallback';
    end

    cycle_results = struct( ...
        'Cycle', {}, ...
        'Session_time_hr', {}, ...
        'Window_start_s', {}, ...
        'Window_end_s', {}, ...
        'Summary_start_s', {}, ...
        'Summary_end_s', {}, ...
        'Fs_Hz', {}, ...
        'Smooth_filter_sec', {}, ...
        'Smooth_filter_order', {}, ...
        'Fit_mode', {}, ...
        'Mean_405_middle3min', {}, ...
        'Mean_465_middle3min', {}, ...
        'Mean_fitted405_middle3min', {}, ...
        'DFF_summary', {}, ...
        'DFF_unsmoothed_summary', {}, ...
        'Is_zscore_reference_cycle', {}, ...
        'Smoothing_fitted405_applied', {}, ...
        'Smoothing_dFF_applied', {}, ...
        'Smoothing_note', {}, ...
        'Fit_slope', {}, ...
        'Fit_intercept', {}, ...
        'Fit_R2', {}, ...
        'Fit_RMSE', {}, ...
        'Fitted405_min_abs_middle3min', {}, ...
        'Fitted405_nonpositive_points_middle3min', {});

    zscore_reference_cycle_means = [];
    zscore_reference_cycles_used = [];

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

        if strcmpi(fit_mode_used, 'reference_cycles')
            p = reference_fit;
        else
            p = polyfit(f405(valid_fit), f465(valid_fit), 1);
        end

        fitted405_raw = polyval(p, f405);
        [fit_r2, fit_rmse] = compute_fit_qc(f465, fitted405_raw, valid_fit);

        %% ===== Smooth fitted 405 control =====
        [fitted405, smoothing_fitted405_applied, fitted405_smooth_note] = ...
            smooth_trace_zero_phase(fitted405_raw, smooth_filter_order);

        %% ===== Calculate dF/F =====
        dff_percent_unsmoothed = (f465 - fitted405) ./ fitted405 * 100;

        %% ===== Smooth dF/F =====
        [dff_percent, smoothing_dff_applied, dff_smooth_note] = ...
            smooth_trace_zero_phase(dff_percent_unsmoothed, smooth_filter_order);

        %% ===== Middle 3 min summary only =====
        summary_idx = tc >= summary_start_sec & tc < summary_end_sec;

        dff_summary = dff_percent(summary_idx);
        dff_unsmoothed_summary = dff_percent_unsmoothed(summary_idx);
        f405_summary = f405(summary_idx);
        f465_summary = f465(summary_idx);
        fitted405_summary = fitted405(summary_idx);

        valid_summary = isfinite(dff_summary) & isfinite(dff_unsmoothed_summary);

        if sum(valid_summary) < 10
            fprintf('Skipping summary cycle %d channel %s: not enough valid dF/F.\n', cyc, label);
            continue
        end

        dff_summary = dff_summary(valid_summary);
        dff_unsmoothed_summary = dff_unsmoothed_summary(valid_summary);
        f405_summary = f405_summary(valid_summary);
        f465_summary = f465_summary(valid_summary);
        fitted405_summary = fitted405_summary(valid_summary);

        fitted405_min_abs = min(abs(fitted405_summary));
        fitted405_nonpositive_points = sum(fitted405_summary <= 0);

        if fitted405_nonpositive_points > 0
            warning('Cycle %d channel %s has %d nonpositive fitted405 points in the summary window.', ...
                cyc, label, fitted405_nonpositive_points);
        end

        mean_dff_for_reference = mean(dff_summary, 'omitnan');
        is_zscore_reference_cycle = ismember(cyc, zscore_reference_cycles);

        if is_zscore_reference_cycle
            zscore_reference_cycle_means = [zscore_reference_cycle_means; mean_dff_for_reference];
            zscore_reference_cycles_used = [zscore_reference_cycles_used; cyc];
        end

        smoothing_note = sprintf('fitted405:%s;dFF:%s', fitted405_smooth_note, dff_smooth_note);

        cycle_results(end+1).Cycle = cyc;
        cycle_results(end).Session_time_hr = t_start / 3600;
        cycle_results(end).Window_start_s = t_start;
        cycle_results(end).Window_end_s = t_end;
        cycle_results(end).Summary_start_s = summary_start_sec;
        cycle_results(end).Summary_end_s = summary_end_sec;
        cycle_results(end).Fs_Hz = fs;
        cycle_results(end).Smooth_filter_sec = smooth_filter_sec;
        cycle_results(end).Smooth_filter_order = smooth_filter_order;
        cycle_results(end).Fit_mode = string(fit_mode_used);
        cycle_results(end).Mean_405_middle3min = mean(f405_summary, 'omitnan');
        cycle_results(end).Mean_465_middle3min = mean(f465_summary, 'omitnan');
        cycle_results(end).Mean_fitted405_middle3min = mean(fitted405_summary, 'omitnan');
        cycle_results(end).DFF_summary = dff_summary;
        cycle_results(end).DFF_unsmoothed_summary = dff_unsmoothed_summary;
        cycle_results(end).Is_zscore_reference_cycle = is_zscore_reference_cycle;
        cycle_results(end).Smoothing_fitted405_applied = smoothing_fitted405_applied;
        cycle_results(end).Smoothing_dFF_applied = smoothing_dff_applied;
        cycle_results(end).Smoothing_note = string(smoothing_note);
        cycle_results(end).Fit_slope = p(1);
        cycle_results(end).Fit_intercept = p(2);
        cycle_results(end).Fit_R2 = fit_r2;
        cycle_results(end).Fit_RMSE = fit_rmse;
        cycle_results(end).Fitted405_min_abs_middle3min = fitted405_min_abs;
        cycle_results(end).Fitted405_nonpositive_points_middle3min = fitted405_nonpositive_points;

    end

    analyzed_cycles = numel(cycle_results);
    if analyzed_cycles < n_cycles
        warning('Channel %s analyzed %d of %d requested cycles.', label, analyzed_cycles, n_cycles);
    end

    missing_reference_cycles = setdiff(zscore_reference_cycles(:), zscore_reference_cycles_used(:));
    if ~isempty(missing_reference_cycles)
        warning('Channel %s z-score reference is missing requested cycle(s): %s.', ...
            label, format_cycle_list(missing_reference_cycles));
    end

    channel_summary = table();

    zscore_mean_dff = mean(zscore_reference_cycle_means, 'omitnan');
    zscore_std_dff = std(zscore_reference_cycle_means, 'omitnan');

    if isempty(zscore_reference_cycle_means) || ~isfinite(zscore_std_dff) || zscore_std_dff == 0
        zscore_mean_dff = NaN;
        zscore_std_dff = NaN;
        warning('Channel %s has no usable first-hour reference dF/F variability for z-scoring.', label);
    end

    requested_reference_label = ['cycles_' format_cycle_list(zscore_reference_cycles) '_middle3min_cycle_mean'];
    used_reference_label = ['cycles_' format_cycle_list(zscore_reference_cycles_used) '_middle3min_cycle_mean'];
    zscore_reference_method = 'cycle_mean_middle3min';

    for result_idx = 1:numel(cycle_results)

        dff_summary = cycle_results(result_idx).DFF_summary;
        dff_unsmoothed_summary = cycle_results(result_idx).DFF_unsmoothed_summary;

        %% ===== Summary parameters =====
        mean_dff   = mean(dff_summary, 'omitnan');
        median_dff = median(dff_summary, 'omitnan');
        max_dff    = max(dff_summary);
        min_dff    = min(dff_summary);
        std_dff    = std(dff_summary, 'omitnan');

        max_unsmoothed_dff = max(dff_unsmoothed_summary);
        min_unsmoothed_dff = min(dff_unsmoothed_summary);
        std_unsmoothed_dff = std(dff_unsmoothed_summary, 'omitnan');

        z_dff_summary = (dff_summary - zscore_mean_dff) ./ zscore_std_dff;

        mean_z_dff   = mean(z_dff_summary, 'omitnan');
        median_z_dff = median(z_dff_summary, 'omitnan');
        max_z_dff    = max(z_dff_summary);
        min_z_dff    = min(z_dff_summary);
        std_z_dff    = std(z_dff_summary, 'omitnan');

        temp_summary = table( ...
            string(label), ...
            cycle_results(result_idx).Cycle, ...
            cycle_results(result_idx).Session_time_hr, ...
            cycle_results(result_idx).Window_start_s, ...
            cycle_results(result_idx).Window_end_s, ...
            cycle_results(result_idx).Summary_start_s, ...
            cycle_results(result_idx).Summary_end_s, ...
            cycle_results(result_idx).Fs_Hz, ...
            cycle_results(result_idx).Smooth_filter_sec, ...
            cycle_results(result_idx).Smooth_filter_order, ...
            cycle_results(result_idx).Fit_mode, ...
            cycle_results(result_idx).Mean_405_middle3min, ...
            cycle_results(result_idx).Mean_465_middle3min, ...
            cycle_results(result_idx).Mean_fitted405_middle3min, ...
            mean_dff, median_dff, max_dff, min_dff, std_dff, ...
            max_unsmoothed_dff, min_unsmoothed_dff, std_unsmoothed_dff, ...
            mean_z_dff, median_z_dff, max_z_dff, min_z_dff, std_z_dff, ...
            zscore_mean_dff, zscore_std_dff, ...
            string(requested_reference_label), string(used_reference_label), ...
            string(zscore_reference_method), ...
            cycle_results(result_idx).Is_zscore_reference_cycle, ...
            cycle_results(result_idx).Smoothing_fitted405_applied, ...
            cycle_results(result_idx).Smoothing_dFF_applied, ...
            cycle_results(result_idx).Smoothing_note, ...
            cycle_results(result_idx).Fit_slope, ...
            cycle_results(result_idx).Fit_intercept, ...
            cycle_results(result_idx).Fit_R2, ...
            cycle_results(result_idx).Fit_RMSE, ...
            cycle_results(result_idx).Fitted405_min_abs_middle3min, ...
            cycle_results(result_idx).Fitted405_nonpositive_points_middle3min, ...
            'VariableNames', { ...
            'Channel','Cycle','Session_time_hr', ...
            'Window_start_s','Window_end_s', ...
            'Summary_start_s','Summary_end_s', ...
            'Fs_Hz','Smooth_filter_sec','Smooth_filter_order','Fit_mode', ...
            'Mean_405_middle3min','Mean_465_middle3min','Mean_fitted405_middle3min', ...
            'Mean_dFF_percent_middle3min','Median_dFF_percent_middle3min', ...
            'Max_dFF_percent_middle3min','Min_dFF_percent_middle3min', ...
            'STD_dFF_percent_middle3min', ...
            'Max_unsmoothed_dFF_percent_middle3min','Min_unsmoothed_dFF_percent_middle3min', ...
            'STD_unsmoothed_dFF_percent_middle3min', ...
            'Mean_z_dFF_middle3min','Median_z_dFF_middle3min', ...
            'Max_z_dFF_middle3min','Min_z_dFF_middle3min', ...
            'STD_z_dFF_middle3min', ...
            'Zscore_reference_mean_dFF_percent','Zscore_reference_STD_dFF_percent', ...
            'Zscore_reference_requested_label','Zscore_reference_used_label', ...
            'Zscore_reference_method','Is_zscore_reference_cycle', ...
            'Smoothing_fitted405_applied','Smoothing_dFF_applied','Smoothing_note', ...
            'Fit_slope','Fit_intercept','Fit_R2','Fit_RMSE', ...
            'Fitted405_min_abs_middle3min','Fitted405_nonpositive_points_middle3min'});

        channel_summary = [channel_summary; temp_summary];

    end

    all_summary = [all_summary; channel_summary];

    if ~isempty(channel_summary) && width(channel_summary) > 0
        writetable(channel_summary, excel_file, 'Sheet', ['Summary_Channel_' label]);
    end

end

%% ================= EXPORT =================
if isempty(all_summary) || width(all_summary) == 0
    warning('No channel summaries were generated. Output files were not written.');
else
    writetable(all_summary, excel_file, 'Sheet', 'Summary_All');
    save(mat_file, 'all_summary');

    fprintf('\nAnalysis complete.\n');
    fprintf('Summary Excel saved to:\n%s\n', excel_file);
    fprintf('Summary MAT saved to:\n%s\n', mat_file);
end

function [reference_fit, ok] = compute_reference_fit(F405, F465, t, fs, reference_cycles, cycle_period_sec, start_sec, end_sec)
    fit405 = [];
    fit465 = [];

    for cyc_idx = 1:numel(reference_cycles)
        cyc = reference_cycles(cyc_idx);
        t_start = (cyc - 1) * cycle_period_sec;
        ref_idx = t >= (t_start + start_sec) & t < (t_start + end_sec);

        if sum(ref_idx) < fs * 10
            continue
        end

        fit405 = [fit405; F405(ref_idx)];
        fit465 = [fit465; F465(ref_idx)];
    end

    valid_fit = isfinite(fit405) & isfinite(fit465);
    ok = sum(valid_fit) >= 10;

    if ok
        reference_fit = polyfit(fit405(valid_fit), fit465(valid_fit), 1);
    else
        reference_fit = [NaN NaN];
    end
end

function [fit_r2, fit_rmse] = compute_fit_qc(signal, fitted, valid_idx)
    valid_idx = valid_idx & isfinite(signal) & isfinite(fitted);

    if sum(valid_idx) < 3
        fit_r2 = NaN;
        fit_rmse = NaN;
        return
    end

    residual = signal(valid_idx) - fitted(valid_idx);
    fit_rmse = sqrt(mean(residual.^2, 'omitnan'));

    signal_centered = signal(valid_idx) - mean(signal(valid_idx), 'omitnan');
    ss_total = sum(signal_centered.^2, 'omitnan');
    ss_residual = sum(residual.^2, 'omitnan');

    if ss_total == 0
        fit_r2 = NaN;
    else
        fit_r2 = 1 - ss_residual / ss_total;
    end
end

function [smoothed_trace, smoothing_applied, smoothing_note] = smooth_trace_zero_phase(trace, filter_order)
    smoothed_trace = double(trace(:));
    smoothing_applied = false;
    smoothing_note = 'skipped';

    if filter_order <= 1
        smoothing_note = 'skipped_filter_order_le_1';
        return
    end

    if numel(smoothed_trace) <= 3 * filter_order
        warning('Skipping smoothing because the trace is too short for filter order %d.', filter_order);
        smoothing_note = 'skipped_short_trace';
        return
    end

    if any(~isfinite(smoothed_trace))
        warning('Skipping smoothing because the trace contains non-finite values.');
        smoothing_note = 'skipped_nonfinite_trace';
        return
    end

    mean_filter = ones(filter_order, 1) / filter_order;
    smoothed_trace = filtfilt(mean_filter, 1, smoothed_trace);
    smoothing_applied = true;
    smoothing_note = 'applied';
end

function cycle_label = format_cycle_list(cycles)
    cycles = cycles(:)';

    if isempty(cycles)
        cycle_label = 'none';
        return
    end

    parts = arrayfun(@num2str, cycles, 'UniformOutput', false);
    cycle_label = strjoin(parts, '_');
end
