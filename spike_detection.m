clear
dataDir = '/Users/asa/Documents/Data/optogenetic_LN_stim/NP1227-Gal4_ACR1 R26A01-LexA_LexAop-mCD8-GFP_PN/2019-06-21';
expName = '2019-06-21_Var_freq_stim__2-hep_+farnesol_10^-2_8s_490_LED_pulse_50p_ND25_ND3_1.mat';
% cd('Z:\Data\recordings\optogenetic_LN_stim\NP1227-Gal4_ACR1 R26A01-LexA_LexAop-mCD8-GFP_PN\2019-06-21')
% rawData = load('2019-06-21_Var_freq_stim__PO_8s_490_LED_pulse_50p_ND25_ND3_1.mat');
cd(dataDir)
rawData = load(expName);
%% Arrange data
data = rawData.data;
sampRate = rawData.sampRate;
spacer_data = rawData.spacer_data;

% Pull out Vm and I for med filter and spike detection
Vm = squeeze(data(:,3,:));
I = squeeze(data(:,1,:));

spacerVm = squeeze(spacer_data(:,3,:));
spacerI = squeeze(spacer_data(:,1,:));


% Set dimensions
VmSize = size(Vm);
spacerVmSize = size(spacerVm);

% Hard coded 100x gain
Vm = (Vm / 100) * 1e3;
spacerVm = (spacerVm / 100) * 1e3;

% Concatenate spacer and experimental trials to make them easier to work
% with. This could also be done by breaking stuff out into a fn.
Vm = cat(1, spacerVm, Vm);
I = cat(1, spacerI, I);

% Filter current trace
bandpassCutoff = [300 1e3];
I = I(:);
I = bandpass_filter(I, bandpassCutoff, sampRate);
I = I / mad(I, 1);

I = reshape(I, VmSize(1) + spacerVmSize(1), VmSize(2));
%% Detect spikes
spd.minProm = 18;
spd.maxWidth = 1e-3 * sampRate;
spd.minWidth = 0.4e-3 * sampRate;
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
    figure
    for iTrial = 1:size(I, 2)
        
        findpeaks(I(:, iTrial) * -1, 'MinPeakProminence', spd.minProm,...
            'MaxPeakWidth', spd.maxWidth,...
            'MinPeakWidth', spd.minWidth,...
            'MinPeakDistance', spd.minDistance,...
            'Annotate','extents');
        hold on
        plot(Vm(:,iTrial))
        axis tight
        title([rawData.matSaveFile, '  - Trial ' num2str(iTrial), ' / ' num2str(VmSize(2))], 'interpreter', 'none')
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