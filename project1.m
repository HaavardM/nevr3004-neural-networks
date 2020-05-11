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

n_bins_angle = 1000; % choose a reasonable number here, try a few
delta_angle = 2*pi / n_bins_angle;

% compute mutual information for each cell, one at a time

zero_mean_rate = 0;

% repeat for HD data
edgesHD = linspace(0, 2*pi, n_bins_angle+1);
mutualInfo = zeros(0, n_cells);
[occupancy,~,angle_inds] = histcounts(headangle,edgesHD); % got an extra output here, what is it?
for i = 1:n_cells
    % get spikes for current cell
    spikes = cellspikes{i};
    if isempty(spikes)
        disp(i);
        continue
    end
    % remove spike times that are outside the range of tracking times
    spikes = spikes(and(spikes >= startTime, spikes <= stopTime));
    
    % bin spike data: from spike times to number of spikes per unit time of tracking data
    % read documentation to see what's happening here!
    edgesT = linspace(startTime,stopTime,numel(trackingtimes)+1);
    binnedSpikes = histcounts(spikes,edgesT);
    
    % try to understand what's happening here! it's a key step
    for iBin = 1:n_bins_angle
        spikesPerAngle(iBin) = sum(binnedSpikes(angle_inds == iBin));
    end
    
    % compute average firing rate for each HD angle
    firing_rate = spikesPerAngle ./ occupancy;
    
    % convert to Hz (spikes per sec) if you haven't already
    firing_rate = 1000*firing_rate / delta_t;
    % this vector is a TUNING CURVE!
    name = {cellnames(i, cellnames(i, :) ~= ' ')};
    if do_plot
        clf;
        plot(rad2deg(edgesHD(1:end-1)), firing_rate);
        print(sprintf('figs/%s', cellnames(i, cellnames(i, :) ~= ' ')), '-depsc');
    end
    % compute probability density (proportion of time spent at each HD angle)
    probability_density = occupancy ./ sum(occupancy);
    
    % compute average firing rate across all HD angles
    mean_rate = mean(firing_rate);
    if or(mean_rate == 0, isnan(mean_rate))
        error("Zero mean rate");
    end
    
    if or(firing_rate == 0, isnan(firing_rate))
        error("Zero firing rate");
    end
    
    if any(probability_density < 0)
        error("Negative probability");
    end
    
    % compute mutual information between firing rate and HD

    % Assume 0log(0) = 0 (in case firing rate is zero)
    % Only sum together fields with non-zero firing rate
    not_zero_inds = (firing_rate ~= 0);
    if ~isempty(not_zero_inds)
        mutualInfo(i) = firing_rate(not_zero_inds) .* log2(firing_rate(not_zero_inds) / mean_rate) * probability_density(not_zero_inds)';
        infoPerSpike(i) = mutualInfo(i) / mean_rate;
    end
end
