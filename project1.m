close all;
[firing_rate, mean_rate, mutual_info, edges_hd, cell_name, P] = part1("Mouse12-120806_awakedata.mat");
%[firing_rate_28, mean_rate_28, mutual_info_28, edges_hd, cell_name] = part1("Mouse28-140313_awakedata.mat");

for i = 1:size(firing_rate, 2)
   clf;
   x = rad2deg(edges_hd(1:end-1));
   y = firing_rate(:, i)';
   plot(x, y);
   name = cell_name(i, cell_name(i, :) ~= ' ');
   title(sprintf("%s, I = %.3f", name, mutual_info(i)));
   xlabel("Head Angle [\circ]");
   ylabel("Firing Rate [Hz]");
   pause;
end


function [firing_rate, mean_rate, mutual_info, edges_hd, cell_name, probability_density] = part1(filename)
    data = load(filename);
    % first and last timestamp
    start_time = data.trackingtimes(1);
    stop_time = data.trackingtimes(end);
    % how much time passes between each timestamp (sampling rate)
    delta_t = data.trackingtimes(2) - data.trackingtimes(1);
    % total number of cells
    n_cells = numel(data.cellspikes);

    n_bins_angle = 360; % choose a reasonable number here, try a few

    edges_t = linspace(start_time,stop_time,numel(data.trackingtimes)+1);
    edges_hd = linspace(0, 2*pi, n_bins_angle+1);
    [occupancy,~,angle_inds] = histcounts(data.headangle,edges_hd); % got an extra output here, what is it?
    
    %preallocate
    mean_rate = zeros(1, n_cells);
    firing_rate = zeros(n_bins_angle, n_cells);
    mutual_info = zeros(1, n_cells);
    spikes_per_angle = zeros(1, n_bins_angle);
    probability_density = zeros(n_bins_angle, n_cells);
    cell_name = data.cellnames;
    
    for i = 1:n_cells
        % get spikes for current cell
        spikes = data.cellspikes{i};
        if isempty(spikes)
            continue
        end
        % remove spike times that are outside the range of tracking times
        spikes = spikes(and(spikes >= start_time, spikes <= stop_time));
        
        % bin spike data: from spike times to number of spikes per unit time of tracking data
        % read documentation to see what's happening here!
        binnedSpikes = histcounts(spikes,edges_t) / (delta_t / 1000);
        
        % try to understand what's happening here! it's a key step
        for iBin = 1:n_bins_angle
            spikes_per_angle(iBin) = sum(binnedSpikes(angle_inds == iBin));
        end
        
        % compute average firing rate for each HD angle
        firing_rate(:, i) = spikes_per_angle ./ occupancy;
        
        % convert to Hz (spikes per sec) if you haven't already
        % this vector is a TUNING CURVE
        firing_rate(:,i) = firing_rate(:,i);
        % compute probability density (proportion of time spent at each HD angle)
        probability_density(:,i) = occupancy ./ sum(occupancy);
        
        % compute average firing rate across all HD angles
        mean_rate(i) = mean(firing_rate(:,i));
        
        % compute mutual information between firing rate and HD

        % Assume 0log(0) = 0 * -inf = 0 (in case firing rate is zero)
        % Only sum together fields with non-zero firing rate
        not_zero_inds = (firing_rate(:, i) ~= 0);
        if ~isempty(not_zero_inds)
            mutual_info(i) = (firing_rate(not_zero_inds, i) .* log2(firing_rate(not_zero_inds, i) / mean_rate(i)))' * probability_density(not_zero_inds, i);
        end
    end
end