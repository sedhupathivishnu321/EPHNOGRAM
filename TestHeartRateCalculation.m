% EPHNOGRAM: A Simultaneous Electrocardiogram and Phonocardiogram Database
% A sample MATLAB script for reading the ECG-PCG data files and extracting the
% R-peaks from the ECG and the S1/S2 peaks of the PCG channels
%
% Dependencies: This script uses some functions from the Open-Source
% Electrophysiological Toolbox (OSET): https://github.com/alphanumericslab/OSET.git
%
% NOTE: This code has not been optimized for any specific dataset and is
% only provided as proof of concept and a starting point to work with the
% EPHNOGRAM database
%
% CITE:
% 1- The EPHNOGRAM on PhysioNet
% 2- A. Kazemnejad, P. Gordany, and R. Sameni. An open-access simultaneous electrocardiogram and phonocardiogram database. bioRxiv, 2021. DOI: https://doi.org/10.1101/2021.05.17.444563
% 3- R. Sameni, The Open-Source Electrophysiological Toolbox (OSET), v 3.14, URL: https://github.com/alphanumericslab/OSET
%
% By: Reza Sameni
% Email: reza.sameni@gmail.com
% May 2021
%

clear;
close all;
clc;

in_folder = '../MATfiles'; % .mat files folder
D = dir([in_folder '/*.mat']); % list of all mat files

% pre-processing parameters
fs = 1000.0; % post-decimation sampling frequency
for m = 1 : length(D) % NOTE: limit the index range to only analyze a subset of the records
    % LOAD DATA
    in_fname = D(m).name;
    in_fname_full = [in_folder '/' in_fname];
    dat = load(in_fname_full); % load the data file
    
    % DECIMATION
    % decimate the ECG and PCG channel to fs
    ECG = decimate(dat.ECG, round(dat.fs/fs));
    PCG = decimate(dat.PCG, round(dat.fs/fs));
    
    % NOTCH FILTERING THE ECG
    fc = 50.0; % powerline frequency
    Qfactor = 45; % Q-factor of the notch filter
    Wo = fc/(fs/2);  BW = Wo/Qfactor; % nothc filter parameters
    [b,a] = iirnotch(Wo, BW); % design the notch filter
    ECG = filtfilt(b, a, ECG); % zero-phase non-causal filtering
    
    % ECG PRE-PROCESSING
    fl_ECG = 10.0; % lower bandpass cut-off frequency for R-peaks
    fh_ECG = 40.0; % upper bandpass cut-off frequency for R-peaks
    ECG_hp = ECG - LPFilter(ECG, fl_ECG/fs);
    ECG_bp = LPFilter(ECG_hp, fh_ECG/fs);
    sigma = 5.0 * std(ECG_bp); % saturate peaks above k-sigma of the ECG
    ECG_saturated = sigma * tanh(ECG_bp / sigma);
    
    if(skew(ECG_saturated) > 0) % find the ECG polarity based on skewness sign
        ecg_polarity = 1;
    else
        ecg_polarity = 0;
    end
    
    % ECG R-PEAK DETECTION
    % First round of R-peak detection (fixed window length peak search)
    ff0 = 1.4; % approximate heart rate in Hz (just a rough estimate; will be fine-tuned by the algorithm)
    peaks0 = PeakDetection(ECG_saturated, ff0/fs, ecg_polarity);
    I0 = find(peaks0);
    ff1 = median(fs ./ diff(I0));
    
    % Second round of R-peak detection (adaptive window length peak search)
    peaks_ECG = PeakDetectionAdaptiveHR(ECG_saturated, ff1, fs, ecg_polarity);
    I_ECG_peaks = find(peaks_ECG);
    RR_intervals_ecg = diff(I_ECG_peaks); % RR-interval time series in samples
    ff_ecg = fs ./ RR_intervals_ecg;
    HR_ecg = 60.* ff_ecg; % ECG-based heart rate
    
    % PCG PRE-PROCESSING
    fl_PCG = 10.0; % lower bandpass cut-off frequency for PCG peaks
    fh_PCG = 100.0; % upper bandpass cut-off frequency for PCG peaks
    PCG_hp = PCG - LPFilter(PCG, fl_PCG/fs);
    PCG_bp = LPFilter(PCG_hp, fh_PCG/fs);
    wlen = round(0.015 * fs); % PCG envelope detector window length
    PCG_envelope = sqrt(filtfilt(ones(1, wlen), wlen, PCG_bp.^2));
    
    % PCG PEAK DETECTION
    search_wlen = round(0.15 * fs);
    % find two dominant local peaks in the PCG envelope between successive
    % ECG R-peaks (gives both S1 and S2)
    num_peaks = 2;
    peaks_S1S2_PCG = IntermediatePeakDetection(PCG_envelope, peaks_ECG, search_wlen, num_peaks, 1, 1);
    % find the first dominant local peak in the PCG envelope between
    % successive ECG R-peaks (gives S1 only)
    num_peaks = 1;
    peaks_S1_PCG = IntermediatePeakDetection(PCG_envelope, peaks_ECG, search_wlen, num_peaks, 1, 2);
    
    I_PCG_S1_peaks = find(peaks_S1_PCG);
    S1S1_intervals_pcg = diff(I_PCG_S1_peaks); % S1S1-interval time series in samples
    ff_pcg = fs ./ S1S1_intervals_pcg;
    HR_pcg_S1 = 60.* ff_pcg; % PCG-based heart rate based on S1S1 time intervals
    
    peaks_S2_PCG = abs(peaks_S1S2_PCG - peaks_S1_PCG); % find the S2 peaks (exclude S1 from the two-peak time sequence)
    I_PCG_S2_peaks = find(peaks_S2_PCG);
    S2S2_intervals_pcg = diff(I_PCG_S2_peaks);  % S2S2-interval time series in samples
    ff_pcg_S2 = fs ./ S2S2_intervals_pcg;
    HR_pcg_S2 = 60.* ff_pcg_S2; % PCG-based heart rate based on S2S2 time intervals
    
    
    % INTERPOLATE THE HEART RATE SEQUENCES
    t = (0 : length(ECG)-1)/fs;
    % Add start and end points to the HR time sequences before interpolation
    tt_ecg = t;
    if(t(I_ECG_peaks(1)) > t(1))
        tt_ecg = cat(2, t(1), tt_ecg(I_ECG_peaks));
        RR_intervals_ecg = cat(2, RR_intervals_ecg(1), RR_intervals_ecg(1), RR_intervals_ecg);
    end
    if(t(I_ECG_peaks(end)) < t(end))
        tt_ecg = cat(2, tt_ecg, t(end));
        RR_intervals_ecg = cat(2, RR_intervals_ecg, RR_intervals_ecg(end));
    end
    tt_pcg_S1 = t;
    if(t(I_PCG_S1_peaks(1)) > t(1))
        tt_pcg_S1 = cat(2, t(1), tt_pcg_S1(I_PCG_S1_peaks));
        S1S1_intervals_pcg = cat(2, S1S1_intervals_pcg(1), S1S1_intervals_pcg(1), S1S1_intervals_pcg);
    end
    if(t(I_PCG_S1_peaks(end)) < t(end))
        tt_pcg_S1 = cat(2, tt_pcg_S1, t(end));
        S1S1_intervals_pcg = cat(2, S1S1_intervals_pcg, S1S1_intervals_pcg(end));
    end
    tt_pcg_S2 = t;
    if(t(I_PCG_S2_peaks(1)) > t(1))
        tt_pcg_S2 = cat(2, t(1), tt_pcg_S2(I_PCG_S2_peaks));
        S2S2_intervals_pcg = cat(2, S2S2_intervals_pcg(1), S2S2_intervals_pcg(1), S2S2_intervals_pcg);
    end
    if(t(I_PCG_S2_peaks(end)) < t(end))
        tt_pcg_S2 = cat(2, tt_pcg_S2, t(end));
        S2S2_intervals_pcg = cat(2, S2S2_intervals_pcg, S2S2_intervals_pcg(end));
    end
    
    % interpolate the HRs
    RR_interval_ecg_interpolated = interp1(tt_ecg, RR_intervals_ecg, t);
    S1S1_interval_pcg_interpolated = interp1(tt_pcg_S1, S1S1_intervals_pcg, t);
    S2S2_interval_pcg_interpolated = interp1(tt_pcg_S2, S2S2_intervals_pcg, t);
    RR_S1S1_diff = (RR_interval_ecg_interpolated - S1S1_interval_pcg_interpolated) / fs;
    RR_S2S2_diff = (RR_interval_ecg_interpolated - S2S2_interval_pcg_interpolated) / fs;
    beta = 0.05; % saturate differences above this amound of time (s)
    RR_S1S1_diff_sat = beta * tanh(RR_S1S1_diff / beta);
    RR_S2S2_diff_sat = beta * tanh(RR_S2S2_diff / beta);
    
    HR_ecg_interpolated = interp1(t(I_ECG_peaks), [HR_ecg(1), HR_ecg], t);
    HR_pcg_S1_interpolated = interp1(t(I_PCG_S1_peaks), [HR_pcg_S1(1), HR_pcg_S1], t);
    HR_pcg_S2_interpolated = interp1(t(I_PCG_S2_peaks), [HR_pcg_S2(1), HR_pcg_S2], t);
    HR_ecg_pcg_S1_diff = HR_ecg_interpolated - HR_pcg_S1_interpolated;
    HR_ecg_pcg_S2_diff = HR_ecg_interpolated - HR_pcg_S2_interpolated;
    alpha = 10.0; % saturate differences above this number of BPMs
    HR_ecg_pcg_S1_diff_sat = alpha * tanh(HR_ecg_pcg_S1_diff / alpha);
    HR_ecg_pcg_S2_diff_sat = alpha * tanh(HR_ecg_pcg_S2_diff / alpha);
    
    % PLOT THE RESULTS
    figure
    lg = {};
    subplot(511);
    plot(t, ECG); lg = cat(2, lg, {'ECG'});
    hold on
    plot(t, ECG_saturated); lg = cat(2, lg, {'saturated ECG'});
    plot(t, PCG); lg = cat(2, lg, {'PCG'});
    plot(t, PCG_envelope); lg = cat(2, lg, {'PCG envelope'});
    plot(t(I_ECG_peaks), ECG(I_ECG_peaks), 'ro'); lg = cat(2, lg, {'R-peaks'});
    plot(t(I_PCG_S1_peaks), PCG_envelope(I_PCG_S1_peaks), 'gx'); lg = cat(2, lg, {'S1-peaks'});
    plot(t(I_PCG_S2_peaks), PCG_envelope(I_PCG_S2_peaks), 'kx'); lg = cat(2, lg, {'S2-peaks'});
    legend(lg);
    grid
    
    subplot(512);
    lg = {};
    plot(t(I_ECG_peaks), [HR_ecg(1), HR_ecg]); lg = cat(2, lg, {'ECG-based HR'});
    hold on
    plot(t(I_PCG_S1_peaks), [HR_pcg_S1(1), HR_pcg_S1]); lg = cat(2, lg, {'S1-based HR'});
    plot(t(I_PCG_S2_peaks), [HR_pcg_S2(1), HR_pcg_S2]); lg = cat(2, lg, {'S2-based HR'});
    legend(lg);
    ylabel('BPM');
    grid
    
    subplot(513);
    lg = {};
    plot(t, HR_ecg_interpolated); lg = cat(2, lg, {'ECG-based HR'});
    hold on
    plot(t, HR_pcg_S1_interpolated); lg = cat(2, lg, {'S1-based HR'});
    plot(t, HR_pcg_S2_interpolated); lg = cat(2, lg, {'S2-based HR'});
    legend(lg)
    grid
    ylabel('Interpolated HRs (BPM)');
    
    subplot(514);
    lg = {};
    plot(t, HR_ecg_pcg_S1_diff_sat); lg = cat(2, lg, {'RR-S1S1'});
    hold on
    plot(t, HR_ecg_pcg_S2_diff_sat); lg = cat(2, lg, {'RR-S2S2'});
    legend(lg)
    grid
    ylabel('RR-SS HR diff saturated (BPM)');
    
    subplot(515);
    lg = {};
    plot(t, 1000 * RR_S1S1_diff_sat); lg = cat(2, lg, {'RR-S1S1'});
    hold on
    plot(t, 1000 * RR_S2S2_diff_sat); lg = cat(2, lg, {'RR-S1S1'});
    legend(lg)
    grid
    ylabel('RR-SS diff saturated (ms)');
    
    % NOTE: Add a keyboard hit or pause here to check the results!
end
