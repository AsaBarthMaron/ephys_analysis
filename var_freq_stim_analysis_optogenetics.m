close all
clear
% dataDir = '/Users/asa/Documents/Data/optogenetic_LN_stim/NP1227-Gal4_ACR1 R26A01-LexA_LexAop-mCD8-GFP_PN/2019-06-21';
dataDir = 'Z:\Data\recordings\optogenetic_LN_stim\NP1227-Gal4_UAS_CsChrimson R26A01-LexA_LexAop-mCD8-GFP_PN\2019-07-09';
expName = '2019-07-09_Var_freq_stim__2-hep_10^-2_490_LED_pulse_same_waveform_30p_ND25_ND3_1.mat';
expName(end-3:end) = [];
cd(fullfile(dataDir, 'analyzed'));
load([expName '_analyzed.mat']);
pulseType = {'10Hz', '2Hz', '0.5Hz'};
%% Load second block & concatenate
% secondExpName = '2019-06-17_Var_freq_stim__2-hep_10^-1_8s_490_LED_pulse_50p_ND25_ND3_reversed_order_1.mat';
% secondExpName(end-3:end) = [];
% secondData = load([secondExpName '_analyzed.mat']);
% 
% expName = [expName '_' secondExpName];
% 
% psth = cat(2, psth, reverse_mat_order(secondData.psth(:,2:end)));
% spacerPsth = cat(2, spacerPsth, reverse_mat_order(secondData.spacerPsth(:,2:end)));
% spacerSpikeInds = cat(2, spacerSpikeInds, reverse_cell_order(secondData.spacerSpikeInds(2:end)));
% spacerVmFilt = cat(2, spacerVmFilt, reverse_mat_order(secondData.spacerVmFilt(:, 2:end)));
% spikeInds = cat(2, spikeInds, reverse_cell_order(secondData.spikeInds(2:end)));
% VmFilt = cat(2, VmFilt, reverse_mat_order(secondData.VmFilt(:, 2:end)));
% 
% % psth = cat(2, psth, secondData.psth(:,2:end));
% % spacerPsth = cat(2, spacerPsth,secondData.spacerPsth(:,2:end));
% % spacerSpikeInds = cat(2, spacerSpikeInds, secondData.spacerSpikeInds(2:end));
% % spacerVmFilt = cat(2, spacerVmFilt, secondData.spacerVmFilt(:, 2:end));
% % spikeInds = cat(2, spikeInds, secondData.spikeInds(2:end));
% % VmFilt = cat(2, VmFilt, secondData.VmFilt(:, 2:end));
% 
% expSize = size(VmFilt);
% spacerSize = size(spacerVmFilt);
% nTrials = expSize(2);
% 
%%
for iTrialType = 2:2:6
    figure
    subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.01], [0.08 0.14], [0.06 0.06]);

    colorOrder = {'k', 'b'};

    %% Raster
    fig = subplot(3, 1, 1);
    ax1 = gca;

    hold on
    % Plot rasters
    for iTrial = iTrialType:6:nTrials
        yPlotInd = iTrial / 6;

        % Plot LED OFF trial rasters
        sp = [spacerSpikeInds{iTrial}; spikeInds{iTrial} + 7e3];
        fig = quickRaster(sp, yPlotInd, 1, colorOrder{1})

        % Plot LED ON trial rasters
        sp = [spacerSpikeInds{iTrial+1}; spikeInds{iTrial+1} + 7e3];
        fig = quickRaster(sp, yPlotInd + ((nTrials-1)/6), 1, colorOrder{2})

        ax1.YDir = 'reverse';
        axis tight
    end
    % Plot spacer / exp divider line
    plot([spacerSize(1) + timeUnits, spacerSize(1) + timeUnits],...
         [0 nTrials/3], 'color', [216, 82, 24] / 255, 'linewidth', 2);

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
    for i = 2:7
        meanPsth(:,i) = mean(p(:,i:6:end), 2);
    end
    meanPsth = meanPsth * timeUnits;

    % Plot psth and adjust axis settings
    plot(meanPsth(:,iTrialType+1:6:end), colorOrder{2}, 'linewidth', 2);
    hold on
    plot(meanPsth(:,iTrialType:6:end), colorOrder{1}, 'linewidth', 2);
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
    yLimits =[-60 0];
 
    fig = subplot(3, 1, 3);
    ax3 = gca;

    % Concatenate spacer and exp
    vf = cat(1, spacerVmFilt, NaN(2 * timeUnits, nTrials), VmFilt);

    % Plot filtered voltage and axis settings
    plot(mean(vf(:,iTrialType+1:6:end), 2), colorOrder{2}, 'linewidth', 2);
    hold on
    plot(mean(vf(:,iTrialType:6:end), 2), colorOrder{1}, 'linewidth', 2);
    axis tight
    xticks([0:timeUnits:(spacerSize(1) + expSize(1) + 2e3)]);
    xticklabels([0:(spacerSize(1) / timeUnits), 0:(expSize(1) / timeUnits)]);
    ax3.TickDir = 'out';
    ax3.FontSize = 14;
    ax3.Box = 'off';
    title('filtered Vm')
    ylim(yLimits)

    % Plot LED and odor
    plot([8.5e3 16.5e3], [yLimits(2) yLimits(2)] - 1, 'linewidth', 8)
    os = downsample(rawData.odorSignal(:,(iTrialType/2)), 10);
    os(os == 0) = NaN;
    os = cat(1, NaN(7e3, 1),  os);
    plot(os * (yLimits(2) -3), 'k', 'linewidth', 10)
    plot([spacerSize(1) + timeUnits, spacerSize(1) + timeUnits],...
         yLimits,  'color', [216, 82, 24] / 255,  'linewidth', 2);

    % Set figure position 
    set(gcf, 'Position', [0,0,1400,600]);
    linkaxes([ax1, ax2, ax3], 'x')
    
    % Save figure
    print([expName '_', pulseType{iTrialType/2}, '.png'], '-dpng', '-r0')
    fig = gcf;
    savefig(gcf, [expName '_', pulseType{iTrialType/2}, '.fig'])
end
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