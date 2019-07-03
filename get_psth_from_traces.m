function psth = get_psth_from_traces(folderPath, dataFilesPath)
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
raster = NaN(trialDuration(1) * sampRate(1), maxReps, nPulseTypes, nBlocks);

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
    Vm = Vm * 1;
    Vm = reshape(Vm, VmSize(1), VmSize(2) * VmSize(3) * VmSize(4));
%     disp('Starting to filter now...')
%     tic
    VmFilt = medfilt1(Vm, 0.08 * sampRate(iBlock), 'truncate');
%     toc
    % VmThres =
    % for i = 1:size(Vm, 2)
    %     [~, locs] = findpeaks(VmThresh(:, i));
    
    checkFiltTrials = 0;
    if checkFiltTrials
        figure
        for i = 1:size(VmFilt, 2)
            plot(Vm(:,i))
            hold on
            plot(VmFilt(:,i))
            hold off
            pause
        end
    end
    
    VmThresh = Vm - VmFilt;
%     VmThresh(VmThresh < 0) = 0;
%     VmThresh = bsxfun(@rdivide, VmThresh, max(VmThresh));
%     VmThresh(VmThresh < 0.3) = 0;
    VmThresh(VmThresh < 7) = 0;
    
    checkThreshTrials = 0;
    if checkThreshTrials
        figure
        for i = 1:size(VmThresh, 2)
            title(['Block ' num2str(iBlock) ', Trial ' num2str(i)])
            plot(VmThresh(:,i))
            pause
        end
    end
    
    filterRaw = 0;
    if filterRaw
        VmFiltRaw = medfilt1(Vm, 0.0050 * sampRate(iBlock), 'truncate');
        checkFiltRaw = 1;
        if checkFiltRaw
            figure
            for i = 1:size(Vm, 2)
                plot(Vm(:,i))
                hold on
                plot(VmFiltRaw(:,i)) 
                hold off
                pause
            end
        end
        Vm = VmFiltRaw;
    end
    
            
    tmpRaster = zeros(size(Vm));
    
    
    for i = 1:size(VmThresh, 2)
%         VmThresh(:,i) = VmThresh(:,i) ./ max(VmThresh(:,i));
        [~, locs] = findpeaks(VmThresh(:, i), 'MinPeakProminence',12, 'MinPeakWidth', sampRate(iBlock) * 0.0007, 'Annotate','extents');
%         [~, locs] = findpeaks(VmThresh(:, i), 'MinPeakProminence',30, 'MinPeakWidth', sampRate(iBlock) * 0.0004, 'MaxPeakWidth', sampRate(iBlock) * 0.007, 'MinPeakDistance', sampRate(iBlock) * 0.001,'Annotate','extents');
%         [~, locs] = findpeaks(VmThresh(:, i), 'MinPeakHeight',0.3, 'MinPeakWidth', sampRate(iBlock) * 0.0004, 'MaxPeakWidth', sampRate(iBlock) * 0.007, 'MinPeakDistance', sampRate(iBlock) * 0.008,'Annotate','extents');
        tmpRaster(locs, i) = 1;
    end
    
    checkPeaks = 1;
    if checkPeaks
        figure
        for i = 1:size(VmFilt, 2)
%             VmThresh(:,i) = VmThresh(:,i) ./ max(VmThresh(:,i));
            title(['Block ' num2str(iBlock) ', Trial ' num2str(i)])
            findpeaks(VmThresh(:, i), 'MinPeakProminence',12, 'MinPeakWidth', sampRate(iBlock) * 0.0007, 'Annotate','extents');
%             findpeaks(VmThresh(:, i), 'MinPeakHeight',0.3, 'MinPeakWidth', sampRate(iBlock) * 0.0004, 'MaxPeakWidth', sampRate(iBlock) * 0.007, 'MinPeakDistance', sampRate(iBlock) * 0.008,'Annotate','extents');
            pause
        end
    end
    
    tmpRaster = reshape(tmpRaster, VmSize(1), VmSize(2), VmSize(3),  VmSize(4));
    tmpPsth = squeeze(mean(tmpRaster, 2));
    % tmpPsth = tmpRaster * 100;
%     binSize = 0.1 * sampRate(iBlock); 
    binSize = 0.2 * sampRate(iBlock); 

    
    for i = 1:VmSize(3)
        for j = 1:VmSize(4)
            tmpPsth(:,i, j) = quickPSTH(tmpPsth(:, i, j), binSize, 'method', 'hanning');
        end
    end
%     for i = 1:VmSize(3)
%         for j = 1:VmSize(4)
%             tmpPsth(:,i, j) = conv(tmpPsth(:, i, j), ones(binSize, 1), 'same');
%         end
%     end
    VmFilt = medfilt1(Vm, 0.04 * sampRate(iBlock), 'truncate');
    VmFilt = reshape(VmFilt, VmSize(1), VmSize(2), VmSize(3),  VmSize(4));
    meanVm(:, :, iBlock) = squeeze(nanmean(VmFilt,2));
    
    raster(:,maxReps:-1:maxReps-(nReps(iBlock)-1),:, iBlock) = tmpRaster;
    psth(:,:,iBlock) = tmpPsth;
    clearvars -except data dataFiles folderPath raster psth olfCh nOlfCh nReps ...
               trialDuration sampRate nPulseTypes maxReps meanVm
end
