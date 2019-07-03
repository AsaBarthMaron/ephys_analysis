function plot_single(varargin)
% Description goes here

%% Parse function inputs
p = inputParser;
p.KeepUnmatched = false;
p.StructExpand = false;
p.CaseSensitive = false;

p.addRequired('data');
p.addRequired('exp');
p.addRequired('sampRate');

p.addOptional('doSave', 0);
p.addOptional('h',0);
p.addParameter('iTrial', '?');
p.addParameter('stim', 0);
p.addOptional('ylims', 0);

p.parse(varargin{:});
data = p.Results.data;
exp = p.Results.exp;
sampRate = p.Results.sampRate;

doSave = p.Results.doSave;
h = p.Results.h;
iTrial = p.Results.iTrial;
stim = p.Results.stim;
ylims = p.Results.ylims;

%% Create figure handle if one is not passed
if ~isa(h,'function_handle')
    h = @(m,n,p) subtightplot (m, n, p, [0.03 0.03], [0.07 0.08], [0.05 0.02]);
end

%% Get amplifier mode and channel info
[mode, units, ~] = get_200B_mode('data', data);
[~, name2num] = get_channel_identities;
vmCh = name2num.ai('10xVm');
currentCh = name2num.ai('Current');
%% Scale the data
scaledOutput = scale_200B_data(data);
switch mode
    case 'I-clamp'
        fixedOutput = (squeeze(data(:,currentCh,:)) / 100) * 1e3; % Fixed I output is set to 100mV/pA
    case 'V-clamp'
        fixedOutput = (squeeze(data(:,vmCh,:)) / 10) * 1e3; % Fixed V output is set to 10x
end
%% Scale stim
if length(ylims)<=1
    ylims(1) = min(scaledOutput(:));
    ylims(2) = max(scaledOutput(:));
    disp(num2str(ylims(2)))
end
if length(stim) > 1
    stim = (stim - min(stim)) / max(stim);
    stim(stim < 0.5) = NaN;
    stim = (stim * ylims(2)) + 1;
end
ylims(2) = ylims(2) + 2;

%% Plot
h(4,1,1)
trialLength = size(scaledOutput, 1)/sampRate;

ax1 = h(5,1,1:4); 
plot((1/sampRate):(1/sampRate):trialLength,...
     scaledOutput(:,1), 'linewidth', 1.2)
switch mode
    case 'I-clamp'
%         ylabel('Membrane Voltage (mV)')
        ylabel('Vm (mV)')
    case 'V-clamp'
        ylabel('Current (pA)')
end
set(ax1, 'box', 'off', 'xcolor', 'k', 'FontSize', 20);
% set(ax1, 'box', 'off', 'xcolor', 'k', 'FontSize', 40);
xlabel('Seconds')
title({[exp.date, ', ' exp.lineName], ...
       [exp.name, ', ' 'Trial #' num2str(iTrial)]}, ...
       'interpreter', 'none')
% axis([0 trialLength ylims(1) ylims(2)])
axis([0 trialLength -70 10])

if length(stim) > 1 
    hold on
    plot((1/sampRate):(1/sampRate):trialLength,...
     stim, 'k', 'LineWidth', 5)
end

% ax2 = h(5,1,5); 
% plot((1/sampRate):(1/sampRate):trialLength,...
%      fixedOutput(:,1), 'color', [0.38 0.38 0.38], 'linewidth', 1.1)
% switch mode
%     case 'I-clamp'
%         ylabel('(pA)')
%     case 'V-clamp'
%         ylabel('(mV)')
% end
% xlabel('Seconds')
set(gcf, 'color', 'w')
% set(ax2, 'box', 'off', 'FontSize', 20);
% lims = axis;
% 
% linkaxes([ax1 ax2], 'x')
% axis([0 trialLength lims(3) lims(4)])

%% Save
if doSave
    if ~isdir([exp.saveDir 'single_trials'])
        mkdir([exp.saveDir filesep 'single_trials']);
    end
    
end
    
    
    
    
    