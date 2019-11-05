% This was created based off the 'var_freq_stim_analysis_optogenetics' but
% for data with only one stimulus waveform. LED trials are assumed to be
% interleaved starting on the second trial.
close all
clear
% saveDir = '~/Documents/Data/optogenetic_LN_stim/2019-07-15_meta_analysis/NP1227-Gal4_ACR1 R26A01-LexA_LexAop-mCD8-GFP_PN/2019-06-17_1';
dataDir = '/Users/asa/Documents/Data/optogenetic_LN_stim/R78F09-Gal4_ACR1 R26A01-LexA_LexAop-mCD8-GFP_PN/2019-10-31';
saveDir = dataDir;
cd(fullfile(dataDir, 'analyzed'));

dataFiles = dir();
dataFiles = dataFiles(~[dataFiles.isdir]);
[~, I] = sort([dataFiles.datenum]);
dataFiles = dataFiles(I);
% dataFiles([5]) = [];
% dataFiles = dataFiles(10);

iStart = 1;
    
expName = dataFiles(end).name;
load(expName);
%%
figure
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.01], [0.08 0.14], [0.06 0.06]);

colorOrder = {'k', 'b'};

%% Raster
fig = subplot(3, 1, 1);
ax1 = gca;

hold on
% Plot rasters
for iTrial = iStart:nTrials
    yPlotInd = iTrial;

    % Plot LED OFF trial rasters
    sp = [spacerSpikeInds{iTrial}; spikeInds{iTrial} + 7e3];
    fig = quickRaster(sp, yPlotInd, 1, colorOrder{1})

    ax1.YDir = 'reverse';
    axis tight
end
% Plot spacer / exp divider line
plot([spacerSize(1) + timeUnits, spacerSize(1) + timeUnits],...
     [0.9 nTrials+1.1], 'color', [216, 82, 24] / 255, 'linewidth', 2);

% Axis settings
xticks([0:timeUnits:(spacerSize(1) + expSize(1) + 2e3)]);
xticklabels([0:(spacerSize(1) / timeUnits), 0:(expSize(1) / timeUnits)]);
ax1.YColor = 'none';
ax1.XAxisLocation = 'top';
ax1.TickDir = 'out';
ax1.FontSize = 14;
title(expName, 'interpreter', 'none')

%% Psth
fig = subplot(3, 1, 2);
ax2 = gca;

% Calculate mean psth & prep for plotting
p = cat(1, spacerPsth, NaN(2 * timeUnits, nTrials), psth);
meanPsth = mean(p, 2);
meanPsth = meanPsth * timeUnits;

% Plot psth and adjust axis settings
plot(meanPsth, colorOrder{1}, 'linewidth', 2);
hold on
axis tight
ax2.XColor = 'none';
ax2.TickDir = 'out';
ax2.Box = 'off';
ax2.FontSize = 14;
title('psth')

% Plot spacer / exp divider line
plot([spacerSize(1) + timeUnits, spacerSize(1) + timeUnits],...
     [0 max(meanPsth(:))], 'color', [216, 82, 24] / 255, 'linewidth', 2);

%% Filtered voltage
yLimits =[-60 -20];

fig = subplot(3, 1, 3);
ax3 = gca;

% Concatenate spacer and exp
vf = cat(1, spacerVmFilt, NaN(2 * timeUnits, nTrials), VmFilt);

% Plot filtered voltage and axis settings
baseline = mean(mean(vf(3*timeUnits:5*timeUnits, :)));
plot([1 spacerSize(1)], [baseline, baseline], '--', 'color', [0.44 0.74 1])
% plot(vf, 'color', [0.8 0.8 0.8], 'linewidth', 2);
hold on
plot([7*timeUnits size(vf,1)], [baseline, baseline], '--', 'color', [0.44 0.74 1])
plot(mean(vf, 2), colorOrder{1}, 'linewidth', 2);
axis tight
xticks([0:timeUnits:(spacerSize(1) + expSize(1) + 2e3)]);
xticklabels([0:(spacerSize(1) / timeUnits), 0:(expSize(1) / timeUnits)]);
ax3.TickDir = 'out';
ax3.FontSize = 14;
ax3.Box = 'off';
title('filtered Vm')
ylim(yLimits)

% Plot LED and odor
% plot([8.25e3 10.25e3], [yLimits(2) yLimits(2)] - 1, 'linewidth', 8)
ls = downsample(rawData.ledSignal, 10);
ls(ls == 0) = NaN;
ls(ls > 0) = 1;
ls = cat(1, NaN(7e3, 1),  ls);
plot(ls * (yLimits(2) - 1), 'b', 'linewidth', 8)
os = downsample(rawData.odorSignal, 10);
os(os == 0) = NaN;
os = cat(1, NaN(7e3, 1),  os);
plot(os * (yLimits(2) -3), 'k', 'linewidth', 10)
plot([spacerSize(1) + timeUnits, spacerSize(1) + timeUnits],...
     yLimits,  'color', [216, 82, 24] / 255,  'linewidth', 2);

% Set figure position 
set(gcf, 'Position', [0,0,1400,600]);
linkaxes([ax1, ax2, ax3], 'x')

% Save figure
if ~isdir([saveDir filesep 'png']) || ~isdir([saveDir filesep 'fig'])
    mkdir([saveDir filesep 'png']);
    mkdir([saveDir filesep 'fig']);
end
print([saveDir filesep 'png' filesep expName, '.png'], '-dpng', '-r0')
fig = gcf;
savefig(gcf, [saveDir filesep 'fig' filesep expName, '.fig'])
    
function data = reverse_mat_order(data)
    tmpData(:,:,1) = data(:,1:2:end);
    tmpData(:,:,2) = data(:,2:2:end);

    data(:,2:2:end) = tmpData(:,:,1);
    data(:,1:2:end) = tmpData(:,:,2);
end
function data = reverse_cell_order(data)
    tmpData(:,1) = data(1:2:end);
    tmpData(:,2) = data(2:2:end);

    data(2:2:end) = tmpData(:,1);
    data(1:2:end) = tmpData(:,2);
end