clear
cd('Z:\Data\recordings\optogenetic_LN_stim\NP1227-Gal4_ACR1 R26A01-LexA_LexAop-mCD8-GFP_PN\2019-06-21')
load('2019-06-21_Var_freq_stim__PO_8s_490_LED_pulse_50p_ND25_ND3_1.mat');

% Pull out Vm and I for med filter and spike detection
Vm = squeeze(data(:,3,:));
I = squeeze(data(:,1,:));

spacerVm = squeeze(spacer_data(:,3,:));
spacerI = squeeze(spacer_data(:,3,:));


% Set dimensions
VmSize = size(Vm);
spacerVmSize = size(spacerVm);

% Hard coded 100x gain
Vm = (Vm / 100) * 1e3;
spacerVm = (spacerVm / 100) * 1e3;

% Filter current trace
I = I(:);
I = bandpass_filter(I, [300 1e3], sampRate);
I = I / mad(I, 1);

spacerI = spacerI(:);
spacerI = bandpass_filter(spacerI, [300 1e3], sampRate);
spacerI = spacerI / mad(spacerI, 1);

I = reshape(I, VmSize(1), VmSize(2));
spacerI = reshape(spacerI, VmSize(1), VmSize(2));
I = cat(1, spacerI, I);
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

checkPeaks = 0;
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
        title([matSaveFile, '  - Trial ' num2str(iTrial), ' / ' num2str(VmSize(2))], 'interpreter', 'none')
        pause
        hold off
    end
end