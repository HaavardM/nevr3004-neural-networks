clear all; close all; clc;

figure(1);
plotPCA('Mouse12-120806_awakedata.mat');
figure(2);
plotPCA('Mouse28-140313_awakedata.mat');



function scores = plotPCA(filename)
    load(filename)
    % first and last timestamp
    startTime = trackingtimes(1);
    stopTime = trackingtimes(end);
    % total number of cells
    n_cells = numel(cellspikes);

    edgesT = linspace(startTime, stopTime, numel(trackingtimes)+1);
    selected = 0;
    for i = 1:n_cells
        % get spikes for current cell
        spikes = cellspikes{i};
        if length(spikes) < 1000
            continue
        end
        % remove spike times that are outside the range of tracking times
        spikes = spikes(and(spikes >= startTime, spikes <= stopTime));
        
        % bin spike data: from spike times to number of spikes per unit time of tracking data
        % read documentation to see what's happening here!
        selected = selected + 1;
        selectedSpikes(selected, :) = histcounts(spikes,edgesT);
        selectedCellNames{selected} = cellnames(i, :);
    end

    selectedSpikesSmooth = smoothdata(selectedSpikes, 2, 'gaussian', 10);
    selectedSpikesStandardized = (selectedSpikesSmooth - mean(selectedSpikesSmooth, 2)) ./ std(selectedSpikesSmooth,0, 2);
    [~, score, ~, ~, explained] = pca(selectedSpikesStandardized');
    scores = score(:, 1:2)';
    scatter(score(:, 1), score(:,2), 5, rad2deg(headangle), 'filled');
    colormap(parula(20))
    c = colorbar;
    c.Label.String = "Head angle [\circ]";
    xlabel(sprintf("PC1 (%.2f %%)", explained(1)));
    ylabel(sprintf("PC2 (%.2f %%)", explained(2)));
    title(filename(1:find(filename == '-')-1));
end