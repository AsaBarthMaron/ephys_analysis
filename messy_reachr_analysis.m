% load('Z:\Data\recordings\optogenetic_LN_stim\GMR-24C12-gal4\2016-06-21\2016-06-21_2ms_pulse_500ms_ipi_1trial_combined_spike_trials_removed')
%%
sampRate = 10e3;
gain = 1e3/50;
stimImpulse = [ones((0.05 * sampRate),1)*1; zeros((0.45 * sampRate),1)];
stim = [zeros(5 * sampRate, 1); repmat(stimImpulse, 60, 1); zeros(5 * sampRate, 1)];
%% Remove artifact
artifactLocs = [zeros(21, 1); 1; 0; 1; zeros((0.5 * sampRate) - 24, 1)];
artifactLocs = [zeros(5 * sampRate, 1); repmat(artifactLocs, 60, 1); zeros(5 * sampRate, 1)];
for iTrial = 1:size(data, 3)
    data(boolean(artifactLocs),3,iTrial) = median(data(:,3,iTrial));
end
%% Filter
Hd = butterworth_for_reachr_vclamp;
for iTrial = 1:size(data, 3)
    filtTrials(:, iTrial) = filtfilt(Hd.sosMatrix, Hd.ScaleValues, data(:,3,iTrial));
end
%% Check individual w/ fft
checkFft = 0;
if checkFft
    T = 1/sampRate;
    L = size(data, 1);
    t = (0:L-1)*T;
    h = figure
    for iTrial = 1:size(data, 3)
        subplot(2,1,1)
        plot_odor_trial(h, data(:,3,iTrial) * gain, stim, sampRate)
        %     hold on
        %     plot_odor_trial(h, filtTrials(:,iTrial) * gain, stim, sampRate)
        %     hold off
        axis([0 40 -50 50])
        ylabel('Current (pA)')
        subplot(2,1,2)
        Y = fft(data(:,3,iTrial));
        P2 = abs(Y/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        f = sampRate * (0:(L/2))/L;
        plot(f(1:5000),P1(1:5000))
        title('Single-Sided Amplitude Spectrum of X(t)')
        xlabel('f (Hz)')
        ylabel('|P1(f)|')
        pause
    end
end
%% Check individual
checkIndv = 1;
if checkIndv
    Hd = butterworth_for_reachr_vclamp;
    h = figure
    for iTrial = 1:size(data, 3)
        clf
        %     subplot(2,1,1)
        %     plot_odor_trial(h, data(:,3,iTrial) * gain, stim, sampRate)
        hold on
        plot_odor_trial(h, filtfilt(Hd.sosMatrix, Hd.ScaleValues, data(:,3,iTrial)) * gain, stim, sampRate)
        hold off
%         axis([0 40 -50 50])
        ylabel('Current (pA)')
        pause
    end
end
%%
h = figure
% plot_odor_trial(h, squeeze(mean(data(:,3,1:4), 3)) * gain, stim, sampRate)
hold on
plot_odor_trial(h, squeeze(mean(cgpFiltTrials, 2)) * gain, stim, sampRate)
% axis([0 40 -40 30])
ylabel('Current (pA)')
title('2016-06-23_cell2, UAS-ReaChR x 70A09-Gal4, 2ms_pulse_500ms_ipi_120s_iti_vclamp', 'interpreter', 'none');
%%
