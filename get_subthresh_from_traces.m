function meanVm = get_subthresh_from_traces(folderPath, dataFilesPath)
    load(dataFilesPath);
    dataFiles = flipud(dataFiles);
    nBlocks = length(dataFiles);
    
    data = load([char(folderPath) filesep dataFiles(1).name]);
    data.spacer_data = [];
    data.spacer_daqInfo = [];
    for iBlock = 2:nBlocks
        tmpData = load([char(folderPath) filesep dataFiles(iBlock).name]);
        tmpData.spacer_data = [];
        tmpData.spacer_daqInfo = [];
        data = [data tmpData];
    end
    clear tmpData

%% 
for iBlock = 1:nBlocks
    nOlfCh(iBlock) = data(iBlock).nOlfCh;
    nReps(iBlock) = data(iBlock).nReps;
    olfCh(iBlock) = data(iBlock).olfCh;
    sampRate(iBlock) = data(iBlock).sampRate;
    trialDuration(iBlock) = size(data(iBlock).data,1) / sampRate(iBlock);
end

nPulseTypes = 3;
maxReps = max(nReps);


%% Remap from linear (random/interleaved) trial structure to a sorted structure. 
% iOdor and pulseType are now lists of indices. 
% TODO: Simplify this so that it only takes single channel olfactometer
% data.

for iBlock = 1:nBlocks
    conditions = data(iBlock).conditions;
    ephysData = data(iBlock).data;
    randTrials = data(iBlock).randTrials;

    
    [iOdor, pulseType] = ind2sub([length(olfCh(iBlock)), 3], randTrials);
    VmSize = [size(ephysData,1), nReps(iBlock), nPulseTypes, nOlfCh(iBlock)];
    Vm = NaN(VmSize(1), VmSize(2), VmSize(3), VmSize(4));
    
    iTrialMap = ones(nOlfCh(iBlock) * nPulseTypes,1);
    for iTrial = 1:length(conditions)
        iTrialType = randTrials(iTrial);
        Vm(:,iTrialMap(iTrialType), pulseType(iTrial), iOdor(iTrial)) ...
            = ephysData(:, 3, iTrial);
%         if iBlock == 5 && iTrial >= 26
%             Vm(:,iTrialMap(iTrialType), pulseType(iTrial), iOdor(iTrial)) ...
%                 = NaN(VmSize(1),1);
%         end
        iTrialMap(iTrialType) = iTrialMap(iTrialType) + 1;
    end
    

    % TODO: get true gain from telegraph output.
    Vm = (Vm/100)* 1e3; % Hard coded 100x gain, rescaling to units of mV.
    %% Find spike times
    % VmThresh = Vm(Vm >
%     dSampFactor = 10;
%     Vm = downsample(Vm, dSampFactor);
%     sampRate(iBlock) = sampRate(iBlock)/dSampFactor;
%     VmSize(1) = VmSize(1)/dSampFactor;
%     Vm = Vm * -1;
    Vm = reshape(Vm, VmSize(1), VmSize(2) * VmSize(3) * VmSize(4));

    VmFilt = medfilt1(Vm, 0.04 * sampRate(iBlock), 'truncate');
    VmFilt = reshape(VmFilt, VmSize(1), VmSize(2), VmSize(3),  VmSize(4));
    meanVm(:, :, iBlock) = squeeze(nanmean(VmFilt,2));

end
