% close all
figure
clear
dataDir = '/Users/asa/Documents/Data/optogenetic_LN_stim/R78F09-Gal4_UAS_CsChrimson R26A01-LexA_LexAop-mCD8-GFP_PN/2019-11-18';
expName = '2019-11-18_1s_farnesol_10^-4_1s_100p_490_LED_ ND25_ND3_1.mat';
cd(dataDir)
rawData = load(expName);
data = scale_200B_data(rawData.data);

%%
nTrials = size(data, 3);

for iTrial = 1:nTrials
    x = data(1.75e4:2.9e4,1,iTrial);
%     x = data(:,1,iTrial);
%     x = [zeros(0.5e4,1); data(1.75e4:2.9e4,1,iTrial)+ 47; zeros(0.5e4,1)];
%     x = data(1e4:4e4,1,iTrial);
%     x = data(2e4:3.15e4,1,iTrial);
%     x = data(1.25e4:2.4e4,1,iTrial);
    rng('default')
    sampRate = 10e3;                                % sample frequency (Hz)


    y = fft(x);

    n = length(x);          % number of samples
    f = (0:n-1)*(sampRate/n);     % frequency range
    power = abs(y).^2/n;    % power of the DFT

    spectra(:,iTrial) = power(1:1000);
end


plot(f(1:1000), mean(spectra(1:1000, 2:2:end),2), 'linewidth', 1.5, 'color', 'k')
% plot(f(1:200), mean(spectra(1:200, :),2), 'linewidth', 1.5, 'color', 'k')
hold on
plot(f(1:1000), mean(spectra(1:1000, 3:2:end),2), 'linewidth', 1.5, 'color', 'r')
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