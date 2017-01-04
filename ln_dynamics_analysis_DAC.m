%% Load all data blocks (corresponding to different odor types).
% dataFiles.mat contains one variable, dataFiles, which is first set using the 
% dir() command in the data folder, and then from which non-odor data files
% are removed.
skipLoad = 0;
if ~skipLoad
    clear
    folderPath = 'Z:\Data\recordings\LN_dynamics\GH298-Gal4\2016-07-29';
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
            plot_odor_trial(h, data(iBlock).data(:,3,iTrial) * 10, ...
                            data(iBlock).odorSignal(:, data(iBlock).randTrials(iTrial)), ...
                            data(iBlock).sampRate)
            ylabel('Vm')
            title([data(iBlock).matSaveFile(12:end-6) ' trial ' num2str(iTrial)],...
                  'interpreter', 'none')
            pause
        end
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
    
    Vm = reshape(Vm, VmSize(1), VmSize(2) * VmSize(3) * VmSize(4));
    disp('Starting to filter now...')
    tic
    VmFilt = medfilt1(Vm, 0.08 * sampRate(iBlock), 'truncate');
    toc
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
    VmThresh(VmThresh < 15) = 0;
    
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
        [~, locs] = findpeaks(VmThresh(:, i), 'MinPeakProminence',10, 'MinPeakWidth', sampRate(iBlock) * 0.0007, 'Annotate','extents');
        tmpRaster(locs, i) = 1;
    end
    
    checkPeaks = 0;
    if checkPeaks
%         figure
        for i = 1:size(VmFilt, 2)
            title(['Block ' num2str(iBlock) ', Trial ' num2str(i)])
            findpeaks(VmThresh(:, i), 'MinPeakProminence',10, 'MinPeakWidth', sampRate(iBlock) * 0.0007, 'Annotate','extents');
            pause
        end
    end
    
    tmpRaster = reshape(tmpRaster, VmSize(1), VmSize(2), VmSize(3),  VmSize(4));
    tmpPsth = squeeze(mean(tmpRaster, 2));
    % tmpPsth = tmpRaster * 100;
    binSize = 0.06 * sampRate(iBlock); 

    
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
    raster(:,maxReps:-1:maxReps-(nReps(iBlock)-1),:, iBlock) = tmpRaster;
    psth(:,:,iBlock) = tmpPsth;
    clearvars -except data dataFiles folderPath raster psth olfCh nOlfCh nReps ...
               trialDuration sampRate nPulseTypes maxReps
end
flipdim(raster, 2);
%% Plot things
figure
subplot = @(m,n,p) subtightplot (m, n, p, [0.03 0.01], [0.03 0.03], [0.03 0.03]);
h = subplot(3, length(dataFiles), 1)

for iBlock = 1:length(dataFiles)
    for iPulseType = 1:3
        h = subplot(3, length(dataFiles), ....
                    iBlock + ((iPulseType -1) * length(dataFiles)))
        rasterLocs = 0:1/size(raster, 2):(1-1/size(raster, 2));
        for iTrial = 1:size(raster, 2)
            h = quickRaster(find(raster(:,iTrial,iPulseType,iBlock)), ...
                                        rasterLocs(iTrial), 1/size(raster, 2));
            hold on
        end
        
        
        h = plot(psth(:,iPulseType,iBlock)/max(psth(:)), 'linewidth', 1.5)
        set(gca, 'box', 'off');
        
        odorSignal = data(iBlock).odorSignal(:, iPulseType);
        odorSignal = (odorSignal/max(odorSignal) * 0.05);
        odorSignal(odorSignal == 0 ) = NaN;
        
        plot(odorSignal + 1, 'color', [0.4940    0.1840    0.5560], 'linewidth', 8);
%         axis tight
        axis([ 0 (11 * sampRate(1)) 0 1.1])
        ax = gca;
        ax.YTick = []
        ax.YAxis.Color = 'w';
        ax.FontSize = 15.4;
%         ax.XLabel.Position = ax.XLabel.Position - [-10e4 -10e4 10e4];
    end
end
for iBlock = 1:length(dataFiles)
    h = subplot(3, length(dataFiles), iBlock)
    title(data(iBlock).matSaveFile(12:end-6), 'interpreter', 'none', 'FontSize', 20)
end
% %% Plot things (for Rachel's RO1 - concentration series)
% figure
% subplot = @(m,n,p) subtightplot (m, n, p, [0.03 0.01], [0.02 0.02], [0.03 0.03]);
% h = subplot(2, length(dataFiles), 1)
% 
% for iBlock = 1:length(dataFiles)
%     for iPulseType = 3
%         h = subplot(2, length(dataFiles), ....
%                     iBlock + ((iPulseType -3) * length(dataFiles)))
%         rasterLocs = 0:1/size(raster, 2):(1-1/size(raster, 2));
%         for iTrial = 1:size(raster, 2)
%             h = quickRaster(find(raster(:,iTrial,iPulseType,iBlock)), ...
%                                         rasterLocs(iTrial), 1/size(raster, 2));
%             hold on
%         end
%         
%         
%         h = plot(psth(:,iPulseType,iBlock)/max(psth(:)), 'linewidth', 1)
%         set(gca, 'box', 'off');
%         
%         odorSignal = data(iBlock).odorSignal(:, iPulseType);
%         odorSignal = (odorSignal/max(odorSignal) * 0.05);
%         odorSignal(odorSignal == 0 ) = NaN;
%         
%         plot(odorSignal + 1, 'k', 'linewidth', 3);
%         axis tight
%         ax = gca;
%         ax.YTick = []
%     end
% end
% for iBlock = 1:length(dataFiles)
%     h = subplot(2, length(dataFiles), iBlock)
%     title(data(iBlock).matSaveFile(12:end-6), 'interpreter', 'none')
% end
% % plot(quickPSTH(psth(:,1,1), binSize))
%% Plot things (for Rachel's RO1 - 70A09 line)
% figure
% subplot = @(m,n,p) subtightplot (m, n, p, [0.03 0.01], [0.02 0.02], [0.03 0.03]);
% h = subplot(3, length(rasterM6), 1)
% 
% for iCell = 1:length(rasterM6)
%     for iPulseType = 1:3
%         h = subplot(3, length(rasterM6), ....
%                     iCell + ((iPulseType -1) * length(rasterM6)))
%         rasterLocs = 0:1/size(rasterM6{iCell}, 2):(1-1/size(rasterM6{iCell}, 2));
%         for iTrial = 1:size(rasterM6{iCell}, 2)
%             h = quickRaster(find(rasterM6{iCell}(:,iTrial,iPulseType)), ...
%                                         rasterLocs(iTrial), 1/size(rasterM6{iCell}, 2));
%             hold on
%         end
%         
%         
%         h = area(psthM6(:,iPulseType,iCell)/max(psthM6(:)), 'linewidth', 1)
%         set(gca, 'box', 'off');
%         
%         odorSignal = data(iCell).odorSignal(:, iPulseType);
%         odorSignal = (odorSignal/max(odorSignal) * 0.05);
%         odorSignal(odorSignal == 0 ) = NaN;
%         
%         plot(odorSignal + 1, 'k', 'linewidth', 3);
%         axis tight
%         ax = gca;
%         ax.YTick = []
%     end
% end
% for iCell = 1:length(rasterM6)
%     h = subplot(3, length(rasterM6), iCell)
%     title(data(iCell).matSaveFile(12:end-6), 'interpreter', 'none')
% end
% %% Plot things (for my DAC - 70A09/24C12 examples)
% figure
% subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.01], [0.02 0.05], [0.045 0.03]);
% h = subplot(3, size(psthM4, 3), 1)
% 
% for iCell = 1:size(psthM4, 3)
%     for iPulseType = 1:3
%         h = subplot(3, size(psthM4, 3), ....
%                     iCell + ((iPulseType -1) * size(psthM4, 3)))
% 
%         
%         h = area(psthM4(:,iPulseType,iCell)/max(psthM4(:)) * 100, 'linewidth', 1)
%         set(gca, 'box', 'off');
%         
%         odorSignal = data(iCell).odorSignal(:, iPulseType);
%         odorSignal = (odorSignal/max(odorSignal) * 0.05);
%         odorSignal(odorSignal == 0 ) = NaN;
%         
%         hold on
%         plot(odorSignal + 102, 'k', 'linewidth', 5);
%         axis([0 11e4 0 105])
%         ax = gca;
%         ax.XTick = []
%         ax.YTick = []
%         ax.YAxis.Color = 'w';
%         ax.XAxis.Color = 'w';
%         ax.FontName = 'Calibri';
%         ax.FontSize = 22;
%     end
% end
% for iCell = 1:size(psthM4, 3)
%     h = subplot(3, size(psthM4, 3), iCell)
%     ax = gca;
%     ax.FontName = 'Calibri';
% %     title(data(iCell).matSaveFile(12:end-6), 'interpreter', 'none')
%     
% end
% for iCell = [1 4 7]
%     subplot(3, size(psthM4, 3), iCell)
%     ax = gca;
%     ax.YAxis.Color = 'k';
%     ax.YTick = [0 100];
%     ax.LineWidth = 3;
% end
% subplot(3,3,4)
% % ylabel('spikes/sec')
% subplot(3,3,1)
% title('Cell 1', 'FontSize', 34)
% subplot(3,3,2)
% title('Cell 2', 'FontSize', 34)
% subplot(3,3,3)
% title('Cell 3', 'FontSize', 34)
% hold on
% plot([nan(8e4,1); ones(2e4,1); nan(1e4,1)] * 70, 'k', 'linewidth', 5);