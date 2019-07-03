function [scaledData, units, sciUnits, mode] = scale_200B_data(varargin)
% Takes data from the 200B amplifer and scales it to common units (mV, pA)
% Data can be in two formats. First is passed in as an M x N x P array, with M 
% being trial length (in samples), N being various DAC inputs, and P being 
% trials. In this case scaled input is #3, gain is input #4 and mode is input #6. 
% The other format is M x P, where the data is just scaled output, while
% gain and mode arrays are passed in separately.

p = inputParser;
p.KeepUnmatched = false;
p.StructExpand = false;
p.CaseSensitive = false;

p.addRequired('data');
p.addOptional('gain',0);
p.addOptional('mode',0);

p.parse(varargin{:});
data = p.Results.data;
gain = p.Results.gain;
mode = p.Results.mode;

[~, name2num] = get_channel_identities;
gainCh = name2num.ai('Gain');
modeCh = name2num.ai('Mode');
dataCh = name2num.ai('Scaled output');

if gain == 0 && mode == 0
    gain = squeeze(data(:,gainCh,:));
    mode = squeeze(data(:,modeCh,:));
end

% Look on page 65 on 200B manual. There are caveats to using some of these
% gain values.
gainLookup = [0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500];

for iTrial = 1:size(data, 3)
    % Find appropriate gain factor (scalingFactor)
    gainIndex = round(median(gain(:, :, iTrial)) * 2); % This is a convenient trick that works
    gainIndex(gainIndex == 0) = []; % Probably should remove this later, just put it in real quick because I had turned off the second amplifier
    scalingFactor = gainLookup(gainIndex);
%     scalingFactor = 100;

    % Scale data
    scaledData(:, :, iTrial) = data(:, dataCh, iTrial)./scalingFactor; % Data in volts/nA
end

if size(scaledData, 2) == 1
    scaledData = squeeze(scaledData);
end
scaledData = scaledData .* 1e3; % Get data in mV or pA
    
% Give appropriate units
[mode, units, sciUnits] = get_200B_mode('mode', mode);
end
        