clear
saveDir = '/Users/asa/Documents/Data/optogenetic_LN_stim/R78F09-Gal4_GFP R60F02-LexA_Chrimson_LN/2021-01-07';
expName = '2021-01-07_1s_2-hep_10^-2_Iclamp_fast_1.mat';
cd(fullfile(saveDir));
load(expName);

% Make/Check fig. folders
if ~isdir([saveDir filesep 'png']) || ~isdir([saveDir filesep 'fig'])
    mkdir([saveDir filesep 'png']);
    mkdir([saveDir filesep 'fig']);
end
%% led only
close all
% plot(mean(data(:,3,:),3) * 10,'linewidth',1.5)
% plot(squeeze(data(:,3,2)) * 10,'linewidth',1.2)
plot(squeeze(data(:,3,:)) * 10,'linewidth',1, 'color', [0.7 0.7 0.7])
hold on
plot(mean(data(:,3,:),3) * 10,'linewidth',1.5, 'color', [0.2 0.2 0.2])
% plot(squeeze(data(:,3,1)) * 10,'linewidth',1.5, 'color', [0.2 0.2 0.2])
plot((ledSignal) +7, 'linewidth', 2)
% plot((odorSignal) -12, 'linewidth', 2)
baseline = mean(mean(data(1:sampRate, 3, :))) * 10;
plot([1 size(data, 1)], [baseline, baseline], '--', 'linewidth', 2)
ylim([-50 10])
xlim([1 size(data,1)])
xlabel('time')
ylabel('Vm')
% ylabel('current, pA')
title(expName, 'interpreter', 'none')
set(gca, 'box', 'off', 'fontsize', 20)
set(gcf, 'Position', [0,0,800,500]);

print([saveDir filesep 'png' filesep expName, '_voltage.png'], '-dpng', '-r0')
fig = gcf;
savefig(gcf, [saveDir filesep 'fig' filesep expName, '_voltage.fig'])
%% led, odor, interleaved
close all
plot(squeeze(data(:,3,2:2:end)) * 10,'linewidth',1, 'color', [0.6 0.6 0.6])
hold on
plot(squeeze(data(:,3,3:2:end)) * 10,'linewidth',1, 'color', [0.9 0.6 0.6])
plot(mean(data(:,3,2:2:end),3) * 10,'linewidth',1.5, 'color', [0.2 0.2 0.2])
plot(mean(data(:,3,3:2:end),3) * 10,'linewidth',1.5, 'color', [0.3 0 0])
% plot(squeeze(data(:,3,1)) * 10,'linewidth',1.5, 'color', [0.2 0.2 0.2])
plot((ledSignal) -7, 'linewidth', 2)
plot((odorSignal) -9, 'linewidth', 2)
% plot((odorSignal(:,2)) +8, 'linewidth', 1.2)
baseline = mean(mean(data(1:sampRate, 3, :))) * 10;
plot([1 size(data, 1)], [baseline, baseline], '--')
ylim([-60 10])
xlim([1 size(data,1)])
xlabel('time')
ylabel('Vm')
% ylabel('current, pA')
title(expName, 'interpreter', 'none')
set(gca, 'box', 'off', 'fontsize', 20)
set(gcf, 'Position', [0,0,800,500]);

print([saveDir filesep 'png' filesep expName, '_voltage.png'], '-dpng', '-r0')
fig = gcf;
savefig(gcf, [saveDir filesep 'fig' filesep expName, '_voltage.fig'])
set(gca, 'box', 'off', 'fontsize', 20)
%% current steps 
clf
% plot(mean(data(:,3,:),3) * 10,'linewidth',1.5)
nSteps = 4;
for i = 1:nSteps
    plot(squeeze(data(:,3,i:nSteps:size(data,3))) * 10,'linewidth',1, 'color', [0.8 0.8 0.8])
    hold on
    plot(mean(data(:,3,i:nSteps:size(data,3)), 3) * 10,'linewidth',1.5, 'color', [0.2 0.2 0.2])
end
baseline = mean(mean(data(1:sampRate, 3, :))) * 10;
plot([1 size(data, 1)], [baseline, baseline], '--')
ylim([-100 50])
xlim([1 size(data,1)])
xlabel('time')
ylabel('Vm')
% title('1s 10^-2 2-hep, 5 trials', 'interpreter', 'none')
% title('1s 10^-2 2-hep, example trial', 'interpreter', 'none')
% title('pA steps = [-100, -50, 50, 100]', 'interpreter', 'none')
title('pA steps = [-60, -30, 30, 60]', 'interpreter', 'none')
set(gca, 'box', 'off', 'fontsize', 20)

print([saveDir filesep 'png' filesep expName, '.png'], '-dpng', '-r0')
fig = gcf;
savefig(gcf, [saveDir filesep 'fig' filesep expName, '.fig'])
%% LED, odor, interleaved 
close all
% plot(mean(data(:,3,:),3) * 10,'linewidth',1.5)
for i = 2:7
    plot(squeeze(data(:,3,i:6:size(data,3))) * 10,'linewidth',1, 'color', [0.8 0.8 0.8])
    hold on
    plot(mean(data(:,3,i:6:size(data,3)), 3) * 10,'linewidth',1.5, 'color', [0.2 0.2 0.2])
end
baseline = mean(mean(data(1:sampRate, 3, :))) * 10;
plot([1 size(data, 1)], [baseline, baseline], '--') 
ylim([-50 0])
xlim([1 size(data,1)])
xlabel('time')
ylabel('Vm')
% title('1s 10^-2 2-hep, 5 trials', 'interpreter', 'none')
% title('1s 10^-2 2-hep, example trial', 'interpreter', 'none')
title('1s 10^-2 2-hep, 2s 10% 490nm ND25 ND3', 'interpreter', 'none')
set(gca, 'box', 'off', 'fontsize', 20)
%% odor only
clf
% plot(mean(data(:,3,:),3) * 10,'linewidth',1.5)
plot(squeeze(data(:,3,:)) * 10,'linewidth',1, 'color', [0.8 0.8 0.8])
hold on
plot(squeeze(data(:,3,2)) * 10,'linewidth',1.5, 'color', [0.2 0.2 0.2])
plot((odorSignal) -5, 'linewidth', 1.2)
ylim([-50 0])
xlim([1 size(data,1)])
xlabel('time')
ylabel('Vm')
% title('1s 10^-2 2-hep, 5 trials', 'interpreter', 'none')
% title('1s 10^-2 2-hep, example trial', 'interpreter', 'none')
title('1s 10^-2 farnesol (black - example trial)', 'interpreter', 'none')
set(gca, 'box', 'off', 'fontsize', 20)
% legend({'all trials', 'example trial'})
%% Odor single pulse + LED
figure
% plot(mean(data(:,3, 2:2:end),3) * 10,'linewidth',1.5)
hold on
% plot(mean(data(:,3, 3:2:end),3) * 10,'linewidth',1.5)
plot(squeeze(data(:,3, :)) * 10,'linewidth',1.1)
plot((ledSignal*10) - 10, 'linewidth', 1.2)
plot(odorSignal - 10, 'linewidth', 1.2)
ylim([-60 0])
xlim([1 size(data,1)])
xlabel('time')
ylabel('mean Vm')
% title('1s 10^-2 2-hep,  2s 5% 490nm LED power,  mean of 4 trials', 'interpreter', 'none')
title('1s 10^-2 2-hep,  2s 5% 565nm LED power,  mean of 3 trials', 'interpreter', 'none')
set(gca, 'box', 'off', 'fontsize', 20)
legend({'odor', 'odor+light'})
%% Var waveform freq + LED
wav([2,4,6]) = {'10Hz', '~2Hz', '~0.5Hz'};
for iStim = [2, 4, 6]
    figure
%     iStim = 2;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
    plot(mean(data(:,3, iStim:6:end),3) * 10,'linewidth',1.5)
    hold on
    plot(mean(data(:,3, iStim+1:6:end),3) * 10,'linewidth',1.5)
    plot((ledSignal*10) - 30, 'linewidth', 1.2)
    plot(odorSignal(:,iStim/2) - 30, 'linewidth', 1.2)
    ylim([-110 0])  
    xlim([1 size(data,1)])
    xlabel('time')
    ylabel('mean Vm')
    title([wav{iStim} ' 10^-1 2-hep,  8s 20% 490nm LED power ND25 ND3,  mean of 6 trials'], 'interpreter', 'none')
    % title('~0.5Hz PO,  8s 10% 490nm LED power ND25 ND3,  mean of 5 trials', 'interpreter', 'none')
    set(gca, 'box', 'off', 'fontsize', 20)
    legend({'odor', 'odor+light'})
end
%% Var waveform freq + LED, split trial blocks
dataDir = 'Z:\Data\recordings\optogenetic_LN_stim\NP1227-Gal4_ACR1 R26A01-LexA_LexAop-mCD8-GFP_PN\2019-06-21\';
load([dataDir ...
    '2019-06-21_Var_freq_stim__2-hep_10^-4_8s_490_LED_pulse_50p_ND25_ND3_1.mat'])
second = load([dataDir ...
    '2019-06-21_Var_freq_stim__2-hep_10^-4_8s_490_LED_pulse_50p_ND25_ND3_reversed_order_1.mat']);


tmpData(:,:,1:2:(size(second.data,3)-1)) = second.data(:,:,3:2:end);
tmpData(:,:,2:2:(size(second.data,3)-1)) = second.data(:,:,2:2:end);
data = cat(3, data, [tmpData]);
wav([2,4,6]) = {'10Hz', '~2Hz', '~0.5Hz'};
for iStim = [2, 4, 6]
    figure
%     iStim = 2;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
    plot(mean(data(:,3, iStim:6:end),3) * 10,'linewidth',1.5)
    hold on
    plot(mean(data(:,3, iStim+1:6:end),3) * 10,'linewidth',1.5)
    plot((ledSignal*10) - 30, 'linewidth', 1.2)
    plot(odorSignal(:,iStim/2) - 30, 'linewidth', 1.2)
    ylim([-70 0])
    xlim([1 size(data,1)])
    xlabel('time')
    ylabel('mean Vm')
    title([wav{iStim} ' 10^-2 2-hep,  8s 50% 490nm LED power ND25 ND3,  mean of 8 trials'], 'interpreter', 'none')
%     title([wav{iStim} ' PO,  8s 50% 490nm LED power ND25 ND3,  mean of 4 trials'], 'interpreter', 'none')
    set(gca, 'box', 'off', 'fontsize', 20)
    legend({'odor', 'odor+light'})
end