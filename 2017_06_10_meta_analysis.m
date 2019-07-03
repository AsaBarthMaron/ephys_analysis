analysisFolder = 'Z:\Data\recordings\LN_dynamics\2017-06-10_meta_analysis';
% load([analysisFolder filesep 'dataFiles.mat']);
% load([analysisFolder filesep 'folderPaths.mat']);
load([analysisFolder filesep 'dataFiles_NP1227_12_19_01_03_all.mat']);
load([analysisFolder filesep 'folderPaths_NP1227_12_19_01_03_all.mat']);

%%
nExps = length(dataFiles);
tmpPsth = zeros(110e3, 3, 6);
psth = zeros(110e3, 3, 6, nExps);
for iExp = 5:7%nExps
    tic
    tmpPsth = get_psth_from_traces(folderPaths{iExp}, dataFiles(iExp).name);
    psth(:,:,:,iExp) = tmpPsth;
    toc
    disp(['Psth calc & saved for exp #' num2str(iExp)])
end
save('psth.mat', 'psth')
%% Mean subtract single trials (for subthresh, irrelevant for PSTHs)
for iExp = 1:nExps
    for iBlock = 1:5
        for iPulseType = 1:3
            psth(:,iPulseType, iBlock, iExp) = psth(:,iPulseType, iBlock, iExp) - mean(psth(:,iPulseType, iBlock, iExp));
        end
    end
end
%% PCA
psth = reshape(psth, 110e3 * 3 * 5, nExps);
mat = downsample(psth,10);
% mat = bsxfun(@minus, mat, min(mat));
% mat = bsxfun(@rdivide, mat, max(mat));
mat = zscore(mat');
% [U, S, latent] = pca(mat');
% mat = mat';
[U,S,V] = svd(mat, 'econ');
figure

subplot(2,2,1)
plot((diag(S)/sum(S(:))) * 100, '*', 'MarkerSize', 15)
title('scree plot')
ylabel('% variance explained')
projPc = mat * V;

%%

% close all
R24C12 = [1 28 29 30];
R70A09 = [2 36 37];
NP1227 = [3:15];
NP2426 = [16:19];
R95G08 = [20, 21];
R11H09 = [22];
R70E03 = [23:25];
R71F08 = [26 27];
NG = [31 34 35];
% NG = [31 35];
GH298 = 32;
R46E11 = 33;

figure
% plot(projPc(:,2), projPc(:,3), '*', 'MarkerSize', 15)
hold on
plot(projPc(R24C12,2), projPc(R24C12,3), '*', 'MarkerSize', 15)
plot(projPc(R70A09,2), projPc(R70A09,3), '*', 'MarkerSize', 15)
plot(projPc(NP1227,2), projPc(NP1227,3), '*', 'MarkerSize', 15)
plot(projPc(NP2426,2), projPc(NP2426,3), '*', 'MarkerSize', 15)
plot(projPc(R95G08,2), projPc(R95G08,3), '*', 'MarkerSize', 15)
plot(projPc(R11H09,2), projPc(R11H09,3), '*', 'MarkerSize', 15)
plot(projPc(R70E03,2), projPc(R70E03,3), '*', 'MarkerSize', 15)
plot(projPc(R71F08,2), projPc(R71F08,3), '*k', 'MarkerSize', 15)
plot(projPc(NG,2), projPc(NG,3), '*r', 'MarkerSize', 15)
plot(projPc(GH298,2), projPc(GH298,3), '*m', 'MarkerSize', 15)
plot(projPc(R46E11,2), projPc(R46E11,3), '*g', 'MarkerSize', 15)

legend({'R24C12', 'R70A09', 'NP1227', 'NP2426', 'R95G08', 'R11H09', 'R70E03', 'R71F08', 'NG', 'GH298', 'R46E11'})
%%
vect = reshape(V, 11e3, 3, 5, nExps);
% for iPC = 1:nExps
iPC = 37;
figure
subplot = @(m,n,p) subtightplot (m, n, p, [0.03 0.015], [0.02 0.02], [0.03 0.03]);
h = subplot(3, 5, 1)
pc = vect(:,:,:,iPC);
% pc = reshape(psth(:,iPC), 110e3, 3, 5);
for iBlock = 1:5
    for iPulseType = 1:3
        h = subplot(3, 5, ....
                    iBlock + ((iPulseType -1) * 5))
        h = plot(pc(:,iPulseType,iBlock), 'linewidth', 1.2)
        hold on
        set(gca, 'box', 'off');
        
%         odorSignal = data(iBlock).odorSignal(:, iPulseType);
%         odorSignal = (odorSignal/max(odorSignal) * (max(meanCurrent(:)) + 5));
%         odorSignal(odorSignal == 0 ) = NaN;
%         
%         plot(odorSignal + 1, 'k', 'linewidth', 3);
% %         axis([0 11e4 min(meanCurrent(:)) max(odorSignal)+2])
%         axis([0 11e4 -80  max(odorSignal)+2])
           axis([0 11e3 min(pc(:))  max(pc(:))])

        ax = gca;
%         ax.YTick = []
    hold off
    end
end
for iBlock = 1:5
    h = subplot(3, 5, iBlock)
    title(dataFiles(iPC).name(1:end-14), 'interpreter', 'none')
end

% pause
end
%% Cluster PCs
for iExp = 1:nExps
    names{iExp} = dataFiles(iExp).name(end-18:end-14);
end
pc = U*S;
% pc = pc(iNs, :);
% pc(30:end,:) = [];
% names = names(iNs);
% names(30:end) = [];
[Y, Z, copDist, meanInconsistency] = cluster_goi(pc(:,1:2));
figure
subplot(2,1,1)
% T = cluster(Z,'cutoff',20);
T = cluster(Z,'maxclust',10);
[~,~, iNs] = dendrogram(Z,50);
ax = gca; ax.XTickLabel = T(iNs);
% ax = gca; ax.XTickLabel = iNs;
title('Clustered on input synapses to LNs (PSDs)')
subplot(2,1,2)
bar(pc(iNs,1:2,1))
set(gca, 'xtick', 1:nExps)
set(gca, 'xticklabel', names(iNs))
axis tight
%% plot individual neurons by clustered order
figure
nConc = 5;
% vect = reshape(V, 11e3, 3, 5, nExps);
for iExp = 5:7%nExps
% iPC = 1;
% figure
subplot = @(m,n,p) subtightplot (m, n, p, [0.03 0.015], [0.02 0.02], [0.03 0.03]);
h = subplot(3, nConc, 1)
% pc = V(:,:,:,iPC);
pc = reshape(psth(:,iExp), 110e3, 3, nConc);
pcNew = reshape(psthNew(:,iExp), 110e3, 3, nConc);
for iBlock = 1:nConc
    for iPulseType = 1:3
        h = subplot(3, nConc, ....
                    iBlock + ((iPulseType -1) * nConc))
        h = plot(pc(:,iPulseType,iBlock), 'linewidth', 1.2)
        hold on
        h = plot(pcNew(:,iPulseType,iBlock), 'linewidth', 1.2)
%         set(gca, 'box', 'off');
        
%         odorSignal = data(iBlock).odorSignal(:, iPulseType);
%         odorSignal = (odorSignal/max(odorSignal) * (max(meanCurrent(:)) + 5));
%         odorSignal(odorSignal == 0 ) = NaN;
%         
%         plot(odorSignal + 1, 'k', 'linewidth', 3);
% %         axis([0 11e4 min(meanCurrent(:)) max(odorSignal)+2])
%         axis([0 11e4 -80  max(odorSignal)+2])
%            axis([0 110e3 min(pc(:))  max(pc(:))])

        ax = gca;
%         ax.YTick = []
    hold off
    end
end
for iBlock = 1:nConc
    h = subplot(3, nConc, iBlock)
    title(dataFiles(iExp).name(1:end-14), 'interpreter', 'none')
end

pause
end
%% plot cluster means
figure
% sortedPsth = psth(iNs,:);
clusterMeans = zeros(length(psth), max(T));
for i = 1:max(T)
    clusterMeans(:,i) = mean(psth(:,T == i),2);
end
% vect = reshape(V, 11e3, 3, 5, nExps);
for iClus = 1:max(T)
% iPC = 1;
% figure
subplot = @(m,n,p) subtightplot (m, n, p, [0.03 0.015], [0.02 0.02], [0.03 0.03]);
h = subplot(3, 5, 1)
% pc = V(:,:,:,iPC);
p = reshape(clusterMeans(:,iClus), 110e3, 3, 5);
% pcNew = reshape(psthNew(:,iExp), 110e3, 3, 5);
for iBlock = 1:5
    for iPulseType = 1:3
        h = subplot(3, 5, ....
                    iBlock + ((iPulseType -1) * 5))
        h = plot(p(:,iPulseType,iBlock), 'linewidth', 1.2)
        hold on
%         h = plot(pcNew(:,iPulseType,iBlock), 'linewidth', 1.2)
%         set(gca, 'box', 'off');
        
%         odorSignal = data(iBlock).odorSignal(:, iPulseType);
%         odorSignal = (odorSignal/max(odorSignal) * (max(meanCurrent(:)) + 5));
%         odorSignal(odorSignal == 0 ) = NaN;
%         
%         plot(odorSignal + 1, 'k', 'linewidth', 3);
% %         axis([0 11e4 min(meanCurrent(:)) max(odorSignal)+2])
%         axis([0 11e4 -80  max(odorSignal)+2])
           axis([0 110e3 min(clusterMeans(:))  max(clusterMeans(:))])

        ax = gca;
%         ax.YTick = []
%     hold off
    end
end
for iBlock = 1:5
    h = subplot(3, 5, iBlock)
    title(num2str(iClus), 'interpreter', 'none')
end

pause
end
%% PRETTY plot for multiple neurons, multiple concentrations, single pulse type
% lnInds = [5, 6 10];
lnInds = 1:2;
% lnInds = [2 1]
nCells = length(lnInds);
nConc = 6;
iPulseType = 3;
sampRate = 1e4;
concentrations = {'Paraffin Oil', '10^-^8', '10^-^6', '10^-^4', '10^-^2', '10^-^1'};
dates = {'spikes/s','spikes/s'}; 
impulse = [ones((2 * sampRate),1)*1; zeros((1.58 * sampRate),1)];
odorSignal = [zeros(2 * sampRate, 1);   repmat(impulse, 2, 1); zeros(ceil(1.84 * sampRate), 1)];
odorSignal(odorSignal == 0 ) = NaN;
%
figure
% subplot = @(m,n,p) subtightplot (m, n, p, [0.14 0.03], [0.18 0.18], [0.05 0.01]);
% subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.03], [0.35 0.35], [0.05 0.01]);
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.03], [0.22 0.22], [0.05 0.01]);
% subplot = @(m,n,p) subtightplot (m, n, p, [0.02 0.03], [0.03 0.03], [0.05 0.01]);
h = subplot(nCells, nConc, 1)

for iCell = 1:nCells
    for iConc = 1:nConc
        h = subplot(nCells, nConc, iConc + ((iCell -1) * nConc));
%         rasterLocs = 0:1/size(rasterM2{iCell}, 2):(1-1/size(rasterM2{iCell}, 2));
%         for iTrial = 1:size(rasterM2{iCell}, 2)
%             h = quickRaster(find(rasterM2{iCell}(:,iTrial,iPulseType)), ...
%                                         rasterLocs(iTrial), 1/size(rasterM2{iCell}, 2));
%             hold on
%         end
%         
        
       h =  plot(downsample(psth(1:10e4,iPulseType,iConc,lnInds(iCell)),100)*sampRate(1));
%        plot(metaPsth{iCell}(:,3,iConc));
        set(gca, 'box', 'off', 'xcolor', 'w', 'ycolor', 'w', 'fontsize', 25);
%         odorSignal = data(iCell).odorSignal(:, iPulseType);
%         odorSignal = (odorSignal/max(odorSignal) * 0.05);
%         odorSignal(odorSignal == 0 ) = NaN;
%         
%         plot(odorSignal + 1, 'k', 'linewidth', 3);
        axis([0 10e2 0 60])
        ax = gca;
%         ax.YTick = []
    
    set(gcf, 'color', 'w')
    hold on
    plot(downsample(odorSignal(1:10e4),100) * 60, 'k', 'linewidth', 5);
        ax = gca;
    set(ax,'xcolor', 'k');
    set(ax, 'xticklabel', {'0' '5' '10'});
    end
end
for iConc = 1:nConc
    h = subplot(nCells, nConc, iConc)
    hold on
%     plot(downsample(odorSignal(1:10e4),10) * 55, 'k', 'linewidth', 5);
    title(concentrations{iConc})
    h = subplot(nCells, nConc, iConc + ((nCells - 1) * nConc))
    ax = gca;
    set(ax,'xcolor', 'k');
    set(ax, 'xticklabel', {'0' '5' '10'});
    
end
for iCell = 1:nCells
    iConc = 1;
    h = subplot(nCells, nConc, iConc + ((iCell -1) * nConc));
    set(gca,'ycolor', 'k');
    ylabel(dates{1})
end
%% Investigate {} tuning
conTun = sum(squeeze(psth(2e4:4e4,3,:,:)),1) - sum(squeeze(psth(1:2e4,3,:,:)),1);
conTun = squeeze(conTun)';
% scatter(conTun(:,1), conTun(:,5))
% corr(conTun(:,1), conTun(:,5))
% plot(conTun')
% plot(bsxfun(@minus, conTun, conTun(:,1))')

[B, I] = sort(conTun(:,1));
bar(B)
a = dataFiles(I);
% Idea is to look at morphology based on this PO response