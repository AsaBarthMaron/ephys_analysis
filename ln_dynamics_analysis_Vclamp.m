%% Load all data blocks (corresponding to different odor types).
% dataFiles.mat contains one variable, dataFiles, which is first set using the 
% dir() command in the data folder, and then from which non-odor data files
% are removed.
skipLoad = 0;
if ~skipLoad
    clear
    folderPath = 'Z:\Data\recordings\LN_dynamics\NP1227-gal4\2017-04-10';
    load([folderPath filesep 'dataFiles.mat']);
    dataFiles = flipud(dataFiles);
    nBlocks = length(dataFiles);
    
    data = load([folderPath filesep dataFiles(1).name]);
%     data = rmfield(data, 'spacer_data');
%     data = rmfield(data, 'spacer_daqInfo');
    for iBlock = 2:nBlocks
        tmpData = load([folderPath filesep dataFiles(iBlock).name]);
%         tmpData.data(:,:,7:30) = zeros(11e4,7,24);
        if isfield(tmpData, 'spacer_data')
%             tmpData = rmfield(tmpData, 'spacer_data');
%             tmpData = rmfield(tmpData, 'spacer_daqInfo');
        end
        if isfield(tmpData, 'access_test')
%             tmpData = rmfield(tmpData, 'access_test');
%             tmpData = rmfield(tmpData, 'access_daqInfo');
        end
        data = [data tmpData];
    end
    clear tmpData
end
Hd = butterworth_for_reachr_vclamp;
%% Remove valve artifcats
% This is much easier to do in Vclamp since they are only a single sample

for iBlock = 1:nBlocks
    for iPulseType = 1:3
        tmp = diff(data(iBlock).odorSignal(:,iPulseType)); % Find pulse onsets/offsets
        tmp = tmp * -1; % Invert them
        artifactLocs(:,1) = find(tmp > 0.5); 
        artifactLocs(:,2) = find(tmp < -0.5);
        artifactLocs = artifactLocs + 2; % Always seem to be 2 samples after onset/offset
        data(iBlock).data(artifactLocs(:),3,iPulseType:3:end) = ...
            data(iBlock).data(artifactLocs(:) - 1,3,iPulseType:3:end)
        clear tmp artifactLocs
    end
end
%% Check individual trials
filter = 1;
filtData = data;
checkRawTrials = 1;
if checkRawTrials
    h = @(m,n,p) subtightplot (m, n, p, [0.03 0.03], [0.13 0.12], [0.055 0.015]);
    for iBlock = 2:nBlocks;
        scaledData = scale_200B_data(data(iBlock).data);
        ylims(2) = max(scaledData(:));
        ylims(1) = min(scaledData(:));
        for iTrial = 1:size(data(iBlock).data, 3)
            stim = data(iBlock).odorSignal(:, data(iBlock).randTrials(iTrial));
            if filter
                trial = data(iBlock).data(:,3,iTrial);
                trial = filtfilt(Hd.sosMatrix, Hd.ScaleValues, trial);
                filtData(iBlock).data(:,3,iTrial) = trial; 
            end
            plot_single(filtData(iBlock).data(:,:,iTrial), data(iBlock).exp,...
                        data(iBlock).sampRate,...
                        'iTria', iTrial, 'stim', stim, 'h', h, 'ylims', ylims);
                    
%             plot_odor_trial(h, data(iBlock).data(:,3,iTrial) * 10, ...
%                             data(iBlock).odorSignal(:, data(iBlock).randTrials(iTrial)), ...
%                             data(iBlock).sampRate)
%             ylabel('Vm')
%             title([data(iBlock).matSaveFile(12:end-6) ' trial ' num2str(iTrial)],...
%                   'interpreter', 'none')
%             set(gcf, 'Position', [0, 0, 1920, 600])
            pause
        end
    end
end 
%% Filter
for iBlock = 1:nBlocks
    current(:,:,iBlock) = scale_200B_data(data(iBlock).data);
        for iTrial = 1:size(current, 2)
            current(:, iTrial, iBlock) = filtfilt(Hd.sosMatrix, Hd.ScaleValues, ...
                current(:,iTrial,iBlock));
        end
end
tmp = current;
current = zeros(11e4, 10, 3, nBlocks);
for iPulseType = 1:3
    current(:,:,iPulseType,:) = tmp(:,iPulseType:3:end,:);
end
meanCurrent = squeeze(mean(current, 2));
%%
figure
subplot = @(m,n,p) subtightplot (m, n, p, [0.03 0.015], [0.02 0.02], [0.03 0.03]);
h = subplot(3, length(dataFiles), 1)

for iBlock = 1:length(dataFiles)
    for iPulseType = 1:3
        h = subplot(3, length(dataFiles), ....
                    iBlock + ((iPulseType -1) * length(dataFiles)))
        h = plot(meanCurrent(:,iPulseType,iBlock), 'linewidth', 1.2)
        hold on
        set(gca, 'box', 'off');
        
        odorSignal = data(iBlock).odorSignal(:, iPulseType);
        odorSignal = (odorSignal/max(odorSignal) * (max(meanCurrent(:)) + 5));
        odorSignal(odorSignal == 0 ) = NaN;
        
        plot(odorSignal + 1, 'k', 'linewidth', 3);
%         axis([0 11e4 min(meanCurrent(:)) max(odorSignal)+2])
        axis([0 11e4 -80  max(odorSignal)+2])

        ax = gca;
%         ax.YTick = []
    end
end
for iBlock = 1:length(dataFiles)
    h = subplot(3, length(dataFiles), iBlock)
    title(data(iBlock).matSaveFile(12:end-6), 'interpreter', 'none')
end