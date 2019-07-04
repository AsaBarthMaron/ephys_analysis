figure
% subplot = @(m,n,p) subtightplot (m, n, p, [0.03 0.01], [0.03 0.03], [0.03 0.03]);


fig = subplot(3, 1, 1)
ax = gca;

hold on
for iTrial = 1:nTrials-5
    sp = [spacerSpikeInds{iTrial}; spikeInds{iTrial} + 7e3];
    fig = quickRaster(sp, iTrial, 1, 'k')
    ax.YDir = 'reverse';
    axis tight
end
% 
% fig = subplot(3, 2, 2)
% ax = gca;
% 
% hold on
% for iTrial = 1:nTrials-5
%     fig = quickRaster(spikeInds{iTrial}, iTrial, 1, 'k')
%     ax.YDir = 'reverse';
%     axis tight
% end