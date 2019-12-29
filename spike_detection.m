close all
clear
dataDir = '/Users/asa/Documents/Data/optogenetic_LN_stim/R24C12-Gal4_ACR1 R26A01-LexA_LexAop-mCD8-GFP_PN/2019-12-11';
expName = '2019-12-11_2Hz_2-hep_10^-2_2s_50p_490_LED_ND25_2.mat';
% cd('Z:\Data\recordings\optogenetic_LN_stim\NP1227-Gal4_ACR1 R26A01-LexA_LexAop-mCD8-GFP_PN\2019-06-21')
% rawData = load('2019-06-21_Var_freq_stim__PO_8s_490_LED_pulse_50p_ND25_ND3_1.mat');
cd(dataDir)
rawData = load(expName);
%% Arrange data
% rawData.data(27502,1,:) = rawData.data(27500,1,:);
data = scale_200B_data(rawData.data);
sampRate = rawData.sampRate;
% d = load('/Users/asa/Documents/Data/optogenetic_LN_stim/R78F09-Gal4_UAS_CsChrimson R26A01-LexA_LexAop-mCD8-GFP_PN/2019-11-18/2019-11-18_1s_farnesol_10^-6_2s_100p_490_LED_ ND25_ND3_1.mat');
% % d = load('2019-11-15_1s_490_LED_pulse_100p_ND25_ND3_1.mat')
% rawData.spacer_data = d.spacer_data(:,:,1:size(rawData.data,3));

% Pull out Vm and I for med filter and spike detection
Vm = squeeze(data(:,1,:));                     % 'scale_200B_data' is still set up for dual patch, so just take 1st channel
I = squeeze(rawData.data(:,1,:));

% Check if spacer data was collected, if not then fill with NaNs & do
% bookkeeping
nTrials = size(rawData.data, 3);
try
    spacer_data = scale_200B_data(rawData.spacer_data);
    spacerVm = squeeze(spacer_data(:,1,:));        % 'scale_200B_data' is still set up for dual patch, so just take 1st channel
    spacerI = squeeze(rawData.spacer_data(:,1,:));
    noSpacer = 0;
catch
    spacerLength = 5;   % Hard coding this to match my most common experiments
    rawData.spacer_data = NaN(5 * sampRate, 2, nTrials);
    spacer_data = rawData.spacer_data;
    spacerVm = squeeze(spacer_data(:,2,:));
    spacerI = spacerVm;
    noSpacer = 1;
end

% Set dimensions
VmSize = size(Vm);
spacerVmSize = size(spacerVm);

% Filter current trace
bandpassCutoff = [300 1e3];
if ~noSpacer
    for iTrial = 1:nTrials
        I(:,iTrial) = bandpass_filter(I(:,iTrial), bandpassCutoff, sampRate);
        I(:,iTrial) = I(:,iTrial) / mad(I(:,iTrial), 1);
        spacerI(:,iTrial) = bandpass_filter(spacerI(:,iTrial), bandpassCutoff, sampRate);
        spacerI(:,iTrial) = spacerI(:,iTrial) / mad(spacerI(:,iTrial), 1);
    end
elseif noSpacer
    for iTrial = 1:nTrials
        I(:,iTrial) = bandpass_filter(I(:,iTrial), bandpassCutoff, sampRate);
        I(:,iTrial) = I(:,iTrial) / mad(I(:,iTrial), 1);
    end
end

% Concatenate spacer and experimental trials to make them easier to work
% with. This could also be done by breaking stuff out into a fn.
Vm = cat(1, spacerVm, Vm);
I = cat(1, spacerI, I);
%% Detect spikes
spd.minProm = 16;
% spd.minProm = 12;
spd.maxWidth = 1.5e-3 * sampRate;
spd.minWidth = 0.3e-3 * sampRate;
spd.minDistance = 1.5e-3 * sampRate;


for iTrial = 1:VmSize(2)
    [~, spikeInds{iTrial}] = findpeaks(I(:, iTrial) * -1,...
        'MinPeakProminence', spd.minProm,...
        'MaxPeakWidth', spd.maxWidth,...
        'MinPeakWidth', spd.minWidth,...
        'MinPeakDistance', spd.minDistance,...
        'Annotate','extents');
end

checkPeaks = 1;
if checkPeaks
%     figure
    clf
    for iTrial = 1:size(I, 2)
        
        findpeaks(I(:, iTrial) * -1, 'MinPeakProminence', spd.minProm,...
            'MaxPeakWidth', spd.maxWidth,...
            'MinPeakWidth', spd.minWidth,...
            'MinPeakDistance', spd.minDistance,...
            'Annotate','extents');
        hold on
        plot(Vm(:,iTrial), 'linewidth', 1.1 , 'color', [0 0.5 0])
        axis tight
        title([expName, '  - Trial ' num2str(iTrial), ' / ' num2str(VmSize(2))], 'interpreter', 'none')
        pause
        hold off
    end
end

%% Create PSTH
% 100ms hanning might be too narrow, but for now I'll stick with it
psthVar.binSize = 0.1 * sampRate;
psthVar.method = 'hanning';

psth = zeros(VmSize(1) + spacerVmSize(1), VmSize(2));

for iTrial = 1:size(I, 2)
    psth(spikeInds{iTrial}, iTrial) = 1;
    psth(:, iTrial) = quickPSTH(psth(:,iTrial),...
                                psthVar.binSize,...
                                'method', psthVar.method);
end
%% Median filter Vm
% 40 ms window
medFiltWindow = 0.04 * sampRate;
VmFilt = medfilt1(Vm, medFiltWindow, 'truncate');
%% Organize variables & save

dsFactor = 10;                   % Downsample factor
timeUnits = sampRate / dsFactor; % In units of bins / second

% Break spacer & experimental trials back apart
spacerVmFilt = VmFilt(1:spacerVmSize(1),:);
spacerPsth = psth(1:spacerVmSize(1),:);
VmFilt = VmFilt(spacerVmSize(1)+1:end, :);
psth = psth(spacerVmSize(1)+1:end, :);

% Do the same for spike times
for iTrial = 1:VmSize(2)
    iSpacerSpikes = spikeInds{iTrial} <= spacerVmSize(1);
    spacerSpikeInds{iTrial} = spikeInds{iTrial}(iSpacerSpikes);
    spikeInds{iTrial} = spikeInds{iTrial}(~iSpacerSpikes) - spacerVmSize(1);
    
    % Put spike times into sample downsampled time base
    spacerSpikeInds{iTrial} = spacerSpikeInds{iTrial} / dsFactor;
    spikeInds{iTrial} = spikeInds{iTrial} / dsFactor;
end

% Downsample data
spacerVmFilt = downsample(spacerVmFilt, dsFactor);
spacerPsth = downsample(spacerPsth, dsFactor);
VmFilt = downsample(VmFilt, dsFactor);
psth = downsample(psth, dsFactor);

% Set data size
spacerSize = size(spacerVmFilt);
expSize = size(VmFilt);
nTrials = expSize(2);

% Remove non-metadata from 'rawData'
rawData = rmfield(rawData, {'data', 'spacer_data'});
clearvars -except spacerVmFilt spacerPsth VmFilt psth spikeInds spacerSpikeInds ...
                  bandpassCutoff dsFactor medFiltWindow psthVar ...
                  spd expSize spacerSize timeUnits ...
                  rawData nTrials expName
if ~isdir('analyzed')
    mkdir('analyzed');
end
cd('analyzed')
save([expName(1:end-4) '_analyzed.mat']);