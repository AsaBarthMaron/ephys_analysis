function h = quickRaster(spikeTimes, Y, height, color)

    bottom = .05;
    top    = .95;
    spikeTimes = spikeTimes(:);
    
    xVals = [spikeTimes,spikeTimes,NaN.*spikeTimes]';
    yVals = [bottom.*ones(length(spikeTimes),1),...
             top.*ones(length(spikeTimes),1),...
             bottom.*ones(length(spikeTimes),1)]';
    h = plot(xVals(:),(yVals(:)*height)+Y, 'color', color, 'linewidth',  1.8);