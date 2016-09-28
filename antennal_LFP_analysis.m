%% Load all data blocks (corresponding to different odor types).
% dataFiles.mat contains one variable, dataFiles, which is first set using the 
% dir() command in the data folder, and then from which non-odor data files
% are removed.
skipLoad = 0;
if ~skipLoad
    clear
    folderPath = 'Z:\Data\recordings\antennal_LFP\wild_type\2016-06-08';
    load([folderPath filesep 'dataFiles.mat']);
    dataFiles = flipud(dataFiles);
    nBlocks = length(dataFiles);
    
    data = load([folderPath filesep dataFiles(1).name]);
    for iBlock = 2:nBlocks
        tmpData = load([folderPath filesep dataFiles(iBlock).name]);
        data = [data tmpData];
    end
    clear tmpData
end
%% Check individual trials

checkRawTrials = 0;
if checkRawTrials
    h = figure;
    for iBlock = nBlocks:-1:1
        for iTrial = 1:size(data(iBlock).data, 3)
            plot_odor_trial(h, data(iBlock).data(:,3,iTrial) * 10, data(iBlock).odorSignal(:, data(iBlock).randTrials(iTrial)), data(iBlock).sampRate)
            ylabel('Vm')
            title([data(iBlock).matSaveFile(12:end-6) ' trial ' num2str(iTrial)], 'interpreter', 'none')
            pause
        end
    end
end 

%% Remap from linear (random/interleaved) trial structure to a sorted structure. 
% iOdor and pulseType are now lists of indices. 
% TODO: Simplify this so that it only takes single channel olfactometer
% data.

raster = NaN(11e4, 10, 3, 5);

for iBlock = 1:nBlocks
    conditions = data(iBlock).conditions;
    ephysData = data(iBlock).data;
    nPulseTypes = 3;
    nOlfCh = data(iBlock).nOlfCh;
    nReps = data(iBlock).nReps;
    olfCh = data(iBlock).olfCh;
    randTrials = data(iBlock).randTrials;
    sampRate = data(iBlock).sampRate;
    
    [iOdor, pulseType] = ind2sub([length(olfCh), 3], randTrials);
    lfpSize = [size(ephysData,1), nReps, nPulseTypes, nOlfCh];
    tmpLfp = NaN(lfpSize(1), lfpSize(2), lfpSize(3), lfpSize(4));
    
    iTrialMap = ones(nOlfCh * nPulseTypes,1);
    for iTrial = 1:length(conditions)
        iTrialType = randTrials(iTrial);
        tmpLfp(:,iTrialMap(iTrialType), pulseType(iTrial), iOdor(iTrial)) ...
            = ephysData(:, 3, iTrial);
        iTrialMap(iTrialType) = iTrialMap(iTrialType) + 1;
    end
    
    % TODO: get true gain from telegraph output.
    tmpLfp = (tmpLfp/500)* 1e3; % Hard coded gain, rescaling to units of mV.
    %% Find spike times
    % VmThresh = Vm(Vm >
    tmpLfp = reshape(tmpLfp, lfpSize(1), lfpSize(2) * lfpSize(3) * lfpSize(4));
    lfpOffset = tmpLfp([1:500 end-(500-1):end],:);
    lfpOffset = mean(lfpOffset, 1);
    tmpLfp = bsxfun(@minus, tmpLfp, lfpOffset);
    
    tmpLfp = reshape(tmpLfp, lfpSize(1), lfpSize(2), lfpSize(3),  lfpSize(4));
    lfp(:,:,:,iBlock) = tmpLfp;
    meanLfp(:,:,iBlock) = squeeze(mean(tmpLfp, 2));

    clearvars -except data dataFiles folderPath meanLfp lfp 
end

%% Remove artifact
removeArtifact = 1;
if removeArtifact
    lfpSize = size(lfp);
    lfp = reshape(lfp, lfpSize(1), lfpSize(2) * lfpSize(3) * lfpSize(4));
    filtLfp = medfilt1(lfp, 0.02 * 1e4, 'truncate');
    
    checkFiltTrials = 0;
    if checkFiltTrials
        figure
        for i = 1:size(filtLfp, 2)
            plot(lfp(:,i))
            hold on
            plot(filtLfp(:,i))
            pause
            hold off
        end
    end
    lfp = reshape(lfp, lfpSize(1), lfpSize(2), lfpSize(3), lfpSize(4));
    filtLfp = reshape(filtLfp, lfpSize(1), lfpSize(2), lfpSize(3), lfpSize(4));
    meanLfp = squeeze(mean(filtLfp, 2));
end 
%% Plot things
figure
subplot = @(m,n,p) subtightplot (m, n, p, [0.03 0.01], [0.02 0.02], [0.03 0.03]);
h = subplot(3, length(dataFiles), 1)

for iBlock = 1:length(dataFiles)
    for iPulseType = 1:3
        h = subplot(3, length(dataFiles), iBlock + ((iPulseType -1)* length(dataFiles)))
%         rasterLocs = 0:1/size(raster, 2):(1-1/size(raster, 2));
%         for iTrial = 1:size(raster, 2)
%             h = quickRaster(find(raster(:,iTrial,iPulseType,iBlock)), rasterLocs(iTrial), 1/size(raster, 2));
%             hold on
%         end
%         
        
        h = plot(meanLfp(:,iPulseType,iBlock), 'linewidth', 1)
        set(gca, 'box', 'off');
        hold on
        odorSignal = data(iBlock).odorSignal(:, iPulseType);
        odorSignal = (odorSignal/max(odorSignal) * 0.05);
        odorSignal(odorSignal == 0 ) = NaN;
        
        plot(odorSignal + 1, 'k', 'linewidth', 3);
        ax = gca;
        linkaxes(ax)
        axis([1 size(meanLfp, 1)  -max(abs(meanLfp(:))) max(abs(meanLfp(:)))/2])
%         ax.YTick = []
    end
end
for iBlock = 1:length(dataFiles)
    h = subplot(3, length(dataFiles), iBlock)
    title(data(iBlock).matSaveFile(12:end-6), 'interpreter', 'none')
end
for iBlock = 1:3
    yInd = 1 + (lfpSize(4) * (iBlock-1));
    subplot(3, length(dataFiles), yInd)
    ylabel('mV')
end