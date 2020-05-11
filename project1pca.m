clear all; close all; clc;
filename = "Mouse12-120806_awakedata.mat";
load(filename)

do_plot = false;
% first and last timestamp
startTime = trackingtimes(1);
stopTime = trackingtimes(end);
% how much time passes between each timestamp (sampling rate)
delta_t = trackingtimes(2) - trackingtimes(1);
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
[coefs, score] = pca(selectedSpikesStandardized');
reduced = score(:, 1:2)';
sampleScores = datasample(score(:, 1:2), 100, 'replace', false);
biplot(coefs(:,1:2),'Scores',sampleScores,'VarLabels',selectedCellNames);