
clear all; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mutual information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load("Mouse12-120806_awakedata.mat")

% first and last timestamp
startTime = trackingtimes(1);
stopTime = trackingtimes(end);
% how much time passes between each timestamp (sampling rate)
deltaT = trackingtimes(2) - trackingtimes(1);
% total number of cells
n_cells = numel(cellspikes);

% compute mutual information for each cell, one at a time
for i = 1:1
    
    % get spikes for current cell
    spikes = cellspikes{i};
    % remove spike times that are outside the range of tracking times
    spikes = spikes(and(spikes >= startTime, spikes <= stopTime));
    
    % bin spike data: from spike times to number of spikes per unit time of tracking data
    % read documentation to see what's happening here!
    edgesT = linspace(startTime,stopTime,numel(trackingtimes)+1);
    binnedSpikes = histcounts(spikes,edgesT);
    
    % repeat for HD data
    n_bins_angle = 10; % choose a reasonable number here, try a few
    edgesHD = linspace(0, 2*pi, n_bins_angle+1);
    [occupancy,~,angle_inds] = histcounts(headangle,edgesHD); % got an extra output here, what is it?
    
    % try to understand what's happening here! it's a key step
    for iBin = 1:n_bins_angle
        spikesPerAngle(iBin) = sum(binnedSpikes(angle_inds == iBin));
    end
    
    % compute average firing rate for each HD angle
    firing_rate = spikesPerAngle ./ occupancy;
    
    % convert to Hz (spikes per sec) if you haven't already
    % this vector is a TUNING CURVE!
    firing_rate = firing_rate * 1000;
    
    % compute probability density (proportion of time spent at each HD angle)
    probability_density = occupancy ./ sum(occupancy);
    
    % compute average firing rate across all HD angles
    mean_rate = mean(firing_rate);
    
    % compute mutual information between firing rate and HD
    mutualInfo = (firing_rate .* log2(firing_rate / mean_rate)) * probability_density';% equation 1 from Skaggs et al. 1993

    % store mutual information for this cell
    %MISSING = mutualInfo;
    
end
