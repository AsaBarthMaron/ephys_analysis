function Rinput = get_Rinput(trace, Icmd, sampRate)

d = diff(Icmd);
cmdTime(1) = find(min(d) == d);
cmdTime(2) = find(max(d) == d);

% if cmdTime(2) < cmdTime(1)
%     error('Ioff should be after Ion, are you using a positive current inj?')
% end
if cmdTime(2) < cmdTime(1)
    cmdTime(1) = 0.5 * sampRate;     % Default current start
    cmdTime(2) = 1 * sampRate;       % Default current stop
end

cmdTime = cmdTime ./ (sampRate/1e3); % convert to ms
cmdTime = round(cmdTime);            % round to nearest ms
cmdTime = cmdTime * (sampRate/1e3);  % convert back to samples

Icmd = Icmd - mean(Icmd(cmdTime(2):end));

IcmdMag = mean(Icmd(cmdTime(1):cmdTime(2)));
IcmdMag = IcmdMag / 0.5e-3;          % This is for fixed 100x gain
IcmdMag = round(IcmdMag);            % I use only integer pA values


calcWin = [cmdTime(1) + (300 * (sampRate/1e3)), cmdTime(2)]; % Takes at least 2-300 ms to settle
V = (mean(trace(calcWin(1):calcWin(2))) - median(trace)) * 1e-3;
IcmdMag = IcmdMag * 1e-12;

Rinput = V / IcmdMag;
Rinput = Rinput / 1e6;
end