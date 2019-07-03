function [psth, raster, odorSignal, dataFiles] = get_pn_psth_from_traces(folderPath, dataFilesPath, stimWav, reversedTrials)
load(fullfile(folderPath, dataFilesPath));
dataFiles = flipud(dataFiles);
nBlocks = length(dataFiles);

data = load([char(folderPath) filesep dataFiles(1).name]);
data.spacer_data = [];
data.spacer_daqInfo = [];
data.impulse = [];
for iBlock = 2:nBlocks
    tmpData = load([char(folderPath) filesep dataFiles(iBlock).name]);
    tmpData.spacer_data = [];
    tmpData.spacer_daqInfo = [];
    tmpData.impulse = [];
    data = [data tmpData];
end
clear tmpData

% Delete interleaved LED trials
for iBlock = 1:nBlocks
    delInds = [1, 3:2:size(data(iBlock).data, 3)];
    if reversedTrials && ~mod(iBlock,2)
        delInds = [1, 2:2:size(data(iBlock).data, 3)];
    end
    data(iBlock).data(:,:,delInds) = [];
    data(iBlock).daqInfo(delInds) = [];
    data(iBlock).randTrials(delInds) = [];
    if strcmpi(stimWav, 'square')
        data(iBlock).nReps = data(iBlock).nReps - length(delInds);
        data(iBlock).conditions(delInds) = [];
    end
end
%%
for iBlock = 1:nBlocks
    nOlfCh(iBlock) = data(iBlock).nOlfCh;
    nReps(iBlock) = data(iBlock).nReps;
    olfCh(iBlock) = data(iBlock).olfCh;
    sampRate(iBlock) = data(iBlock).sampRate;
    trialDuration(iBlock) = size(data(iBlock).data,1) / sampRate(iBlock);
end

nPulseTypes = size(data(1).odorSignal,2);
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
    ISize = VmSize;
    I = NaN(ISize(1), ISize(2), ISize(3), ISize(4));
    
    iTrialMap = ones(nOlfCh(iBlock) * nPulseTypes,1);
    for iTrial = 1:length(conditions)
        iTrialType = randTrials(iTrial);
        Vm(:,iTrialMap(iTrialType), pulseType(iTrial), iOdor(iTrial)) ...
            = ephysData(:, 3, iTrial);
        I(:,iTrialMap(iTrialType), pulseType(iTrial), iOdor(iTrial)) ...
            = ephysData(:, 1, iTrial);
        %         if iBlock == 5 && iTrial >= 26
        %             Vm(:,iTrialMap(iTrialType), pulseType(iTrial), iOdor(iTrial)) ...
        %                 = NaN(VmSize(1),1);
        %         end
        iTrialMap(iTrialType) = iTrialMap(iTrialType) + 1;
    end
    
    
    % TODO: get true gain from telegraph output.
    Vm = (Vm/100)* 1e3; % Hard coded 100x gain, rescaling to units of mV.
    %%
    I = I(:);
%     I = lowPassFilter(I, 1e3, sampRate(iBlock));
    I = bandpass_filter(I, [300 1e3], sampRate(iBlock));    % Bandpass filter
    I = I / mad(I, 1);  % Normalize to median average deviation
    I = reshape(I, ISize(1), ISize(2) * ISize(3) * ISize(4));
    tmpRaster = zeros(size(I));
    
    
    for i = 1:size(I, 2)
        [~, locs] = findpeaks(I(:, i)*-1, 'MinPeakProminence', 13, 'MaxPeakWidth', 1e-3 * sampRate(iBlock), 'MinPeakWidth', 0.4e-3 * sampRate(iBlock), 'MinPeakDistance', 1.5e-3 * sampRate(iBlock),  'Annotate','extents');
        tmpRaster(locs, i) = 1;

%         % Find positive * negative peaks
%         [~, locsUs] = findpeaks(zscore(I(:, i)*-1), 'MinPeakHeight', 3, 'MinPeakDistance', 1.5e-3 * sampRate(iBlock),  'Annotate','extents');
%         [~, locsDs] = findpeaks(zscore(I(:, i)), 'MinPeakHeight', 4, 'MinPeakDistance', 1.5e-3 * sampRate(iBlock),  'Annotate','extents');
%         
%         % Write both pos and neg peak locations to tmp raster
%         tmpRaster(locsUs, i) = 1;
%         tmpRaster(locsDs, i) = 1;
%         % Eliminate 
%         delInds = [];
%         for iLoc = 1:length(locsUs)
%             sumWindow = 2.2e-3 * sampRate(iBlock);
%             windowEnd = locsUs(iLoc) + sumWindow;
%             if windowEnd > length(locsUs)
%                 windowEnd = iLoc;
%             end
%             nLobes = sum(tmpRaster(locsUs(iLoc):windowEnd, i));
%             if nLobes == 1
%                 delInds = [delInds, iLoc];
%             elseif nLobes ~= 2
%                 error(['More than 2 lobes detected within ' num2str(sumWindow/10) 'ms'])
%             end
%         end
%         
%         tmpRaster(:,i) = zeros(ISize(1), 1);
%         locsUs(delInds) = [];
%         tmpRaster(locsUs, i) = 1;
%                 
        
    end
    
    checkPeaks = 1;
    if checkPeaks
        figure
%         for i = 1:size(I, 2)
%             zI = zscore(I(:,i));
%             plot(zI, 'linewidth', 1.1)
%             hold on
%             locs = find(tmpRaster(:,i));
%             plot(locs, zI(locs), 'O', 'linewidth', 2)
%             hold off
%             pause
%         end
        for i = 1:size(I, 2)
            %             VmThresh(:,i) = VmThresh(:,i) ./ max(VmThresh(:,i));
%                         findpeaks(I(:, i), 'MinPeakProminence',3, 'MinPeakWidth', sampRate(iBlock) * 0.7e-3, 'Annotate','extents');
            findpeaks(I(:, i)*-1, 'MinPeakProminence', 13, 'MaxPeakWidth', 1e-3 * sampRate(iBlock), 'MinPeakWidth', 0.4e-3 * sampRate(iBlock), 'MinPeakDistance', 1.5e-3 * sampRate(iBlock),  'Annotate','extents');
            hold on
%             findpeaks(zscore(I(:, i)*-1), 'MinPeakHeight', 4, 'MinPeakDistance', 1.5e-3 * sampRate(iBlock),  'Annotate','extents');
            %             findpeaks(VmThresh(:, i), 'MinPeakHeight',0.3, 'MinPeakWidth', sampRate(iBlock) * 0.0004, 'MaxPeakWidth', sampRate(iBlock) * 0.007, 'MinPeakDistance', sampRate(iBlock) * 0.008,'Annotate','extents');
            title(['Block ' num2str(iBlock) ', Trial ' num2str(i)])
            pause
            hold off
        end
    end
    
    tmpRaster = reshape(tmpRaster, VmSize(1), VmSize(2), VmSize(3),  VmSize(4));
    tmpPsth = squeeze(mean(tmpRaster, 2));
    % tmpPsth = tmpRaster * 100;
    binSize = 0.1 * sampRate(iBlock);
    %     binSize = 0.2 * sampRate(iBlock);
    
    
    for i = 1:ISize(3)
        for j = 1:ISize(4)
            tmpPsth(:,i, j) = quickPSTH(tmpPsth(:, i, j), binSize, 'method', 'hanning');
        end
    end
    raster(:,maxReps:-1:maxReps-(nReps(iBlock)-1),:, iBlock) = tmpRaster;
    psth(:,:,iBlock) = tmpPsth;
    odorSignal = data(1).odorSignal;
    clearvars -except data dataFiles folderPath raster psth odorSignal olfCh nOlfCh nReps ...
        trialDuration sampRate nPulseTypes maxReps meanVm
        %% Find spike times
%         % VmThresh = Vm(Vm >
%     %     dSampFactor = 10;
%     %     Vm = downsample(Vm, dSampFactor);
%     %     sampRate(iBlock) = sampRate(iBlock)/dSampFactor;
%     %     VmSize(1) = VmSize(1)/dSampFactor;
%         Vm = Vm * 1;    
%         Vm = Vm(:);
%         Vm = lowPassFilter(Vm, 1e3, sampRate(iBlock));
% 
%         Vm = reshape(Vm, VmSize(1), VmSize(2) * VmSize(3) * VmSize(4));
%     %     disp('Starting to filter now...')
%     %     tic
%         VmFilt = medfilt1(Vm, 0.08 * sampRate(iBlock), 'truncate');
%     %     toc
%         % VmThres =
%         % for i = 1:size(Vm, 2)
%         %     [~, locs] = findpeaks(VmThresh(:, i));
%     
%         checkFiltTrials = 0;
%         if checkFiltTrials
%             figure
%             for i = 1:size(VmFilt, 2)
%                 plot(Vm(:,i))
%                 hold on
%                 plot(VmFilt(:,i))
%                 hold off
%                 pause
%             end
%         end
%     
%         VmThresh = Vm - VmFilt;
%     %     VmThresh(VmThresh < 0) = 0;
%     %     VmThresh = bsxfun(@rdivide, VmThresh, max(VmThresh));
%     %     VmThresh(VmThresh < 0.3) = 0;
%         VmThresh(VmThresh < 7) = 0;
%     
%         checkThreshTrials = 0;
%         if checkThreshTrials
%             figure
%             for i = 1:size(VmThresh, 2)
%                 title(['Block ' num2str(iBlock) ', Trial ' num2str(i)])
%                 plot(VmThresh(:,i))
%                 pause
%             end
%         end
%     
%         filterRaw = 0;
%         if filterRaw
%             VmFiltRaw = medfilt1(Vm, 0.0050 * sampRate(iBlock), 'truncate');
%             checkFiltRaw = 1;
%             if checkFiltRaw
%                 figure
%                 for i = 1:size(Vm, 2)
%                     plot(Vm(:,i))
%                     hold on
%                     plot(VmFiltRaw(:,i))
%                     hold off
%                     pause
%                 end
%             end
%             Vm = VmFiltRaw;
%         end
%     
%     
%         tmpRaster = zeros(size(Vm));
%     
%     
%         for i = 1:size(VmThresh, 2)
%     %         VmThresh(:,i) = VmThresh(:,i) ./ max(VmThresh(:,i));
%             [~, locs] = findpeaks(VmThresh(:, i), 'MinPeakProminence',12, 'MinPeakWidth', sampRate(iBlock) * 0.0007, 'Annotate','extents');
%     %         [~, locs] = findpeaks(VmThresh(:, i), 'MinPeakProminence',30, 'MinPeakWidth', sampRate(iBlock) * 0.0004, 'MaxPeakWidth', sampRate(iBlock) * 0.007, 'MinPeakDistance', sampRate(iBlock) * 0.001,'Annotate','extents');
%     %         [~, locs] = findpeaks(VmThresh(:, i), 'MinPeakHeight',0.3, 'MinPeakWidth', sampRate(iBlock) * 0.0004, 'MaxPeakWidth', sampRate(iBlock) * 0.007, 'MinPeakDistance', sampRate(iBlock) * 0.008,'Annotate','extents');
%             tmpRaster(locs, i) = 1;
%         end
%     
%         checkPeaks = 1;
%         if checkPeaks
%             figure
%             for i = 1:size(VmFilt, 2)
%     %             VmThresh(:,i) = VmThresh(:,i) ./ max(VmThresh(:,i));
%                 title(['Block ' num2str(iBlock) ', Trial ' num2str(i)])
%                 findpeaks(Vm(:,i) - VmFilt(:, i), 'MinPeakProminence',12, 'MinPeakWidth', sampRate(iBlock) * 0.0007, 'Annotate','extents');
%     %             findpeaks(VmThresh(:, i), 'MinPeakHeight',0.3, 'MinPeakWidth', sampRate(iBlock) * 0.0004, 'MaxPeakWidth', sampRate(iBlock) * 0.007, 'MinPeakDistance', sampRate(iBlock) * 0.008,'Annotate','extents');
%                 pause
%             end
%         end
%     
%         tmpRaster = reshape(tmpRaster, VmSize(1), VmSize(2), VmSize(3),  VmSize(4));
%         tmpPsth = squeeze(mean(tmpRaster, 2));
%         % tmpPsth = tmpRaster * 100;
%     %     binSize = 0.1 * sampRate(iBlock);
%         binSize = 0.2 * sampRate(iBlock);
%     
%     
%         for i = 1:VmSize(3)
%             for j = 1:VmSize(4)
%                 tmpPsth(:,i, j) = quickPSTH(tmpPsth(:, i, j), binSize, 'method', 'hanning');
%             end
%         end
%     %     for i = 1:VmSize(3)
%     %         for j = 1:VmSize(4)
%     %             tmpPsth(:,i, j) = conv(tmpPsth(:, i, j), ones(binSize, 1), 'same');
%     %         end
%     %     end
%         VmFilt = medfilt1(Vm, 0.04 * sampRate(iBlock), 'truncate');
%         VmFilt = reshape(VmFilt, VmSize(1), VmSize(2), VmSize(3),  VmSize(4));
%         meanVm(:, :, iBlock) = squeeze(nanmean(VmFilt,2));
%     
%         raster(:,maxReps:-1:maxReps-(nReps(iBlock)-1),:, iBlock) = tmpRaster;
%         psth(:,:,iBlock) = tmpPsth;
%         clearvars -except data dataFiles folderPath raster psth olfCh nOlfCh nReps ...
%                    trialDuration sampRate nPulseTypes maxReps meanVm
end
