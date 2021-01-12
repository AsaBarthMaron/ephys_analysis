close all
clear
dataDir = '/Users/asa/Documents/Data/optogenetic_LN_stim/R78F09-Gal4_GFP R60F02-LexA_Chrimson_LN/2020-11-24';
expName = '2020-11-24_1s_2-hep_10^-2_2s_490_LED_pulse_100p_1.mat'
cd(dataDir)
rawData = load(expName);
data = scale_200B_data(rawData.data);

%%
nTrials = size(data, 3);

for iTrial = 1:nTrials
    x = data(1.75e4:2.9e4,1,iTrial);
%     x = data(1.25e4:2.4e4,1,iTrial);
    rng('default')
    sampRate = 10e3;                                % sample frequency (Hz)


    y = fft(x);

    n = length(x);          % number of samples
    f = (0:n-1)*(sampRate/n);     % frequency range
    power = abs(y).^2/n;    % power of the DFT

    spectra(:,iTrial) = power(1:200);
end


plot(f(1:200), mean(spectra(1:200, 2:2:end),2), 'linewidth', 1.5, 'color', 'k')
% plot(f(1:200), mean(spectra(1:200, :),2), 'linewidth', 1.5, 'color', 'k')
hold on
plot(f(1:200), mean(spectra(1:200, 3:2:end),2), 'linewidth', 1.5, 'color', 'r')
xlabel('Frequency')
ylabel('Power')
ylim([0 2000])
xlim([0 200])
hold on
% title ('Odor (2-heptanone) evoked frequency spectrum: CSDn > Chrimson, lLN2F_b > GFP', 'interpreter', 'none')
title ('PO evoked frequency spectrum: CSDn > Chrimson, lLN2F_b > GFP', 'interpreter', 'none')
legend({'light off', 'light on'})
% print([dataDir filesep 'png' filesep expName, '_spectral_analysis.png'], '-dpng', '-r0')
% fig = gcf;
% savefig(gcf, [dataDir filesep 'fig' filesep expName, '_spectral_analysis.fig'])
%% Figure
% 
% figure
% plot(0:(1/sampRate):1.7, data(1.5e4:3.2e4,1,2), 'linewidth', 1.3, 'color', 'k')
% hold on
% plot(0:(1/sampRate):1.7, data(1.5e4:3.2e4,1,3), 'linewidth', 1.3, 'color', 'r')
% title ('Odor (2-heptanone): CSDn > Chrimson, lLN2F_b > GFP', 'interpreter', 'none')
% legend({'light off', 'light on'})