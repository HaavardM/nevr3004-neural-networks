close all;

%Task 1b


filename_12 = "Mouse12-120806_awakedata.mat";
filename_28 = "Mouse28-140313_awakedata.mat";
%[firing_rate, mean_rate, mutual_info, edges_hd, cell_name, P] = part1("Mouse12-120806_awakedata.mat");
[firing_rate_12, mean_rate_12, mutual_info_12, edges_hd, cell_name_12, P_12] = part1(filename_12);
[firing_rate_28, mean_rate_28, mutual_info_28, edges_hd, cell_name_28, P_28] = part1(filename_28);
cell_name_12 = cellstr(cell_name_12);
cell_name_28 = cellstr(cell_name_28);
brain_area_12 = strings(14,1);
brain_area_12(1:8) = "Thalamus";
brain_area_12(14) = "Hippocampus";

brain_area_28 = strings(14,1);
brain_area_28(1:7) = "Post-Subiculum";
brain_area_28(8:11) = "Thalamus";


%1A
x = rad2deg(edges_hd(1:end-1));
% for i = 1:size(firing_rate_12, 2)
%    clf;
%    name = cell_name_12{i};
%    grp = str2double(name(2));
%    name = "Mouse 12: " + name + sprintf(" (%s)", brain_area_12(grp));
%    plotTuning(x, firing_rate_12(:, i), mutual_info_12(i), name)
%    print(sprintf("report/project1/figs/%s/%s", filename_12, cell_name_12{i}), "-depsc");
%    pause;
% end
% 
% for i = 1:size(firing_rate_28, 2)
%    clf;
%    name = cell_name_28{i};
%    grp = str2double(name(2));
%    name = "Mouse 28: " + name + sprintf(" (%s)", brain_area_28(grp));
%    plotTuning(x, firing_rate_28(:, i), mutual_info_28(i), name)
%    print(sprintf("report/project1/figs/%s/%s", filename_28, cell_name_28{i}), "-depsc");
%    pause;
% end

%1B
C12_T3C10 = find(cell_name_12 == "T3C10");
C12_T4C7 = find(cell_name_12 == "T4C7");
C28_T8C7 = find(cell_name_28 == "T8C7");
C28_T4C14 = find(cell_name_28 == "T4C14");

figure(2);
subplot(2,2,1);
plot(x, firing_rate_12(:, C12_T3C10));
title({"Mouse 12: T3C10", "Thalamus", sprintf("I = %.2f", mutual_info_12(C12_T3C10))});
xlabel("Head angle [\circ]"); ylabel("Firing Rate [Hz]");
xlim([0, 360]);
subplot(2,2,2);
plot(x, firing_rate_12(:, C12_T4C7));
title({"Mouse 12: T4C7", "Thalamus", sprintf("I = %.2f", mutual_info_12(C12_T4C7))});
xlabel("Head angle [\circ]"); ylabel("Firing Rate [Hz]");
xlim([0, 360]);
subplot(2,2,3);
plot(x, firing_rate_28(:, C28_T8C7));
title({"Mouse 28: T8C7", "Thalamus", sprintf("I = %.2f", mutual_info_28(C28_T8C7))});
xlabel("Head angle [\circ]"); ylabel("Firing Rate [Hz]");
xlim([0, 360]);
subplot(2,2,4);
plot(x, firing_rate_28(:, C28_T4C14));
title({"Mouse 28: T4C14", "Post Subiculum", sprintf("I = %.2f", mutual_info_28(C28_T4C14))});
xlabel("Head angle [\circ]"); ylabel("Firing Rate [Hz]");
xlim([0, 360]);
print("report/project1/figs/1b", "-depsc");

%

figure(3);
subplot(1,2,1);
plot(x, P_12);
title({"Mouse 12", "P(X)"});
xlim([0, 360]);
xlabel("Head direction [\circ]");
ylabel("Probability");
subplot(1,2,2);
plot(x, P_28);
title({"Mouse 28", "P(X)"});
xlabel("Head direction [\circ]");
ylabel("Probability");
xlim([0, 360]);
print("report/project1/figs/probs", "-depsc");

figure(4);
C12_T4C10 = find(cell_name_12 == "T4C10");
C28_T3C7 = find(cell_name_28 == "T3C7");
C28_T9C8 = find(cell_name_28 == "T9C8");
subplot(3,1,1);
name = cell_name_12{C12_T4C10};
grp = str2double(name(2));
name = "Mouse 12: " + name + sprintf(" (%s)", brain_area_12(grp));
plotTuning(x, firing_rate_12(:, C12_T4C10), mutual_info_12(C12_T4C10), name);
subplot(3,1,2);
name = cell_name_28{C28_T3C7};
grp = str2double(name(2));
name = "Mouse 28: " + name + sprintf(" (%s)", brain_area_28(grp));
plotTuning(x, firing_rate_28(:, C28_T3C7), mutual_info_28(C28_T3C7), name);
subplot(3,1,3);
name = cell_name_28{C28_T9C8};
grp = str2double(name(2));
name = "Mouse 28: " + name + sprintf(" (%s)", brain_area_28(grp));
plotTuning(x, firing_rate_28(:, C28_T9C8), mutual_info_28(C28_T9C8), name);
print("report/project1/figs/peculiar_tuning", "-depsc");

figure(5);
C28_T4C3 = find(cell_name_28 == "T4C3");
C28_T4C4 = find(cell_name_28 == "T4C4");
C28_T4C11 = find(cell_name_28 == "T4C11");
subplot(3,1,1);
name = cell_name_28{C28_T4C3};
grp = str2double(name(2));
name = "Mouse 28: " + name + sprintf(" (%s)", brain_area_28(grp));
plotTuning(x, firing_rate_28(:, C28_T4C3), mutual_info_28(C28_T4C3), name);
subplot(3,1,2);
name = cell_name_28{C28_T4C4};
grp = str2double(name(2));
name = "Mouse 28: " + name + sprintf(" (%s)", brain_area_28(grp));
plotTuning(x, firing_rate_28(:, C28_T4C4), mutual_info_28(C28_T4C4), name);
subplot(3,1,3);
name = cell_name_28{C28_T4C11};
grp = str2double(name(2));
name = "Mouse 28: " + name + sprintf(" (%s)", brain_area_28(grp));
plotTuning(x, firing_rate_28(:, C28_T4C11), mutual_info_28(C28_T4C11), name);
print("report/project1/figs/weird_mutual", "-depsc");

function plotTuning(x,y, mut_info, cell_name)
   plot(x, y);
   name = cell_name;
   title({name, sprintf("I = %.3f",mut_info)});
   xlabel("Head Angle [\circ]");
   ylabel("Firing Rate [Hz]");
   xlim([0, 360]);
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
    probability_density = occupancy ./ sum(occupancy);
    
    %preallocate
    mean_rate = zeros(1, n_cells);
    firing_rate = zeros(n_bins_angle, n_cells);
    mutual_info = zeros(1, n_cells);
    spikes_per_angle = zeros(1, n_bins_angle);
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
        
        
        % compute average firing rate across all HD angles
        mean_rate(i) = mean(firing_rate(:,i));
        
        % compute mutual information between firing rate and HD

        % Assume 0log(0) = 0 * -inf = 0 (in case firing rate is zero)
        % Only sum together fields with non-zero firing rate
        not_zero_inds = (firing_rate(:, i) ~= 0);
        if ~isempty(not_zero_inds)
            mutual_info(i) = (firing_rate(not_zero_inds, i) .* log2(firing_rate(not_zero_inds, i) / mean_rate(i)))' * probability_density(not_zero_inds)';
        end
    end
    %Mutual information below zero does not make any sense
    mutual_info(mutual_info < 0) = 0;
end