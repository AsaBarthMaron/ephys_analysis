% This was created based off the 'var_freq_stim_analysis_optogenetics' but
% for data with only one stimulus waveform. LED trials are assumed to be
% interleaved starting on the second trial.
close all
clear
% saveDir = '~/Documents/Data/optogenetic_LN_stim/2019-07-15_meta_analysis/NP1227-Gal4_ACR1 R26A01-LexA_LexAop-mCD8-GFP_PN/2019-06-17_1';
dataDir = '/Users/asa/Documents/Data/optogenetic_LN_stim/R78F09-Gal4_GFP R60F02-LexA_Chrimson_LN/2020-12-15';
saveDir = dataDir;
cd(fullfile(dataDir, 'analyzed'));

dataFiles = dir();  
dataFiles = dataFiles(~[dataFiles.isdir]);
[~, I] = sort([dataFiles.datenum]);
dataFiles = dataFiles(I);
% dataFiles([5]) = [];
% dataFiles = dataFiles(10);

iStart = 2;
% for iExp = 1:length(dataFiles)    
    
    expName = dataFiles(end).name;
    load(expName);
%%

isi_control = [];
isi_light = [];

for iTrial = 2:nTrials
    r = rand(1);
    if ~mod(iTrial, 2)
        isi_control = [isi_control; diff(spikeInds{iTrial})];
    elseif mod(iTrial,2)
        isi_light = [isi_light; diff(spikeInds{iTrial})];
    end
end

edges = 0:3:100;
figure
histogram(isi_control, edges, 'normalization', 'probability', 'facecolor', 'k');
hold on
histogram(isi_light, edges, 'normalization', 'probability', 'facecolor', [0.7 0 0 ]);
xlabel('Inter-spike interval (ms)')
ylabel('Probability')
set(gca, 'box', 'off', 'tickdir', 'out', 'fontsize', 20)


%% Save fig
if ~isdir([saveDir filesep 'png']) || ~isdir([saveDir filesep 'fig'])
    mkdir([saveDir filesep 'png']);
    mkdir([saveDir filesep 'fig']);
end
print([saveDir filesep 'png' filesep expName '_ISI_hist' '.png'], '-dpng', '-r0')
fig = gcf;
savefig(gcf, [saveDir filesep 'fig' filesep expName '_ISI_hist' '.fig'])