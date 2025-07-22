% Clear workspace and figures
clear;
close all;

% Paths for the input data that you would like to compare
tpath_ori ='';
tpath_ori_2 ='';

% Parameters
xlimit_1 = -100;
xlimit_2 = 200;
freq = 1200;
binx = (xlimit_1:1:xlimit_2-1); % Define bin range

%% Obtain trace data for the first set

% Change directory to the first path
cd(tpath_ori);
raw_ori = dir('*'); %bring in the data
num_file_ori = numel(raw_ori); % Number of files in the folder

% Initialize cell arrays and matrices to store data
dz_ori = cell(1, num_file_ori);
med_dz_ori = cell(1, num_file_ori);
dz_collect_ori_med = [];

% Load data from each file and apply median filtering
for i = 1:num_file_ori
    dz_ori{i} = load(raw_ori(i).name); % Load data from each text file
    med_dz_ori{i} = medfilt1(dz_ori{i}, freq/10); % Apply median filter
end

% Collect all data into single arrays
for k = 1:num_file_ori
    dz_collect_ori_med = [dz_collect_ori_med, med_dz_ori{k}']; % Concatenate filtered data
end

% Plot histogram for the first set of data
figure('visible', 'off');
h_raw_ind_all = histogram(dz_collect_ori_med, 'BinWidth', 1, 'BinLimits', [xlimit_1, xlimit_2]);
raw_bar_all = h_raw_ind_all.Values';
Total_Time_all = (raw_bar_all / num_file_ori) / freq;

% Plot the normalized time histogram
figure(88124);
bar(binx, Total_Time_all, 1, 'EdgeColor', [0.50, 0.50, 0.50], 'FaceAlpha', 0.5);
hold on;

%% Obtain trace data for the second set

% Change directory to the second path
cd(tpath_ori_2);
raw_ori_2 = dir('*');
num_file_ori_2 = numel(raw_ori_2); % Number of files in the second folder

% Initialize cell arrays and matrices for the second data set
dz_ori_2 = cell(1, num_file_ori_2);
med_dz_ori_2 = cell(1, num_file_ori_2);
dz_collect_ori_med_2 = [];

% Load and filter data for the second set
for i = 1:num_file_ori_2
    dz_ori_2{i} = load(raw_ori_2(i).name); % Load data
    med_dz_ori_2{i} = medfilt1(dz_ori_2{i}, freq/10); % Apply median filter
end

% Collect all data into single arrays
for k = 1:num_file_ori_2
    dz_collect_ori_med_2 = [dz_collect_ori_med_2, med_dz_ori_2{k}']; % Concatenate filtered data
end

% Plot histogram for the second set of data
figure('visible', 'off');
h_raw_ind_all_2 = histogram(dz_collect_ori_med_2, 'BinWidth', 1, 'BinLimits', [xlimit_1, xlimit_2]);
raw_bar_all_2 = h_raw_ind_all_2.Values';
Total_Time_all_2 = (raw_bar_all_2 / num_file_ori_2) / freq;

% Plot the normalized time histogram
figure(88124);
bar(binx, Total_Time_all_2, 1, 'EdgeColor', [0.50, 0.50, 0.50], 'FaceAlpha', 0.5);

%% Random comparing xlimits designated

% Set the number of runs
num_runs = 100;
how_many_random = 40; % Number of random selections

% Get user input for range
to_where_start = input('At what point do you want to start? ');
to_where_end = input('At what point do you want to end? ');

% Convert user input to indices within the histogram
new_x_s = to_where_start - xlimit_1 + 1;
new_x_e = to_where_end - xlimit_1 + 1;

% Initialize arrays to store results
p_val_1 = zeros(num_runs, 1);
p_val_2 = zeros(num_runs, 1);

% Run random selection and calculation for specified number of runs
for run = 1:num_runs
    % Select random indices
    random_indices_all = randperm(num_file_ori, how_many_random);
    random_indices_all_2 = randperm(num_file_ori_2, how_many_random);
    
    % Collect random samples from each data set
    dz_collect_ori_r = [];
    dz_collect_ori_r_2 = [];
    for i = 1:how_many_random
        dz_collect_ori_r = [dz_collect_ori_r, med_dz_ori{random_indices_all(i)}'];
        dz_collect_ori_r_2 = [dz_collect_ori_r_2, med_dz_ori_2{random_indices_all_2(i)}'];
    end

    % Calculate histograms for random samples
    figure('visible', 'off');
    h_raw_ind = histogram(dz_collect_ori_r, 'BinWidth', 1, 'BinLimits', [xlimit_1, xlimit_2]);
    h_raw_ind_2 = histogram(dz_collect_ori_r_2, 'BinWidth', 1, 'BinLimits', [xlimit_1, xlimit_2]);
    
    % Normalize the histogram data
    raw_bar = h_raw_ind.Values';
    Total_Time = (raw_bar / how_many_random) / freq;
    raw_bar_2 = h_raw_ind_2.Values';
    Total_Time_2 = (raw_bar_2 / how_many_random) / freq;
    
    % Calculate area under curve for the specified range
    p_val_1(run) = sum(Total_Time(new_x_s:new_x_e));
    p_val_2(run) = sum(Total_Time_2(new_x_s:new_x_e));
end

% Calculate means of the area under curve
mean_p_val_1 = mean(p_val_1);
mean_p_val_2 = mean(p_val_2);

% Perform a t-test between the two sets of AUC values
[h, p] = ttest2(p_val_1, p_val_2, 'Alpha', 0.01);

% Create data for box plot
data_for_boxplot = [p_val_1, p_val_2];

% Create box plot for AUC comparisons
figure;
boxplot(data_for_boxplot, 'Labels', {'+', '-'});
ylabel('Values');
title(['AUC of ', num2str(to_where_start), ' to ', num2str(to_where_end), ' nm']);

% Overlay individual data points
hold on;
scatter(ones(size(p_val_1)), p_val_1, 'r', 'filled');
scatter(ones(size(p_val_2)) * 2, p_val_2, 'b', 'filled');

% Determine significance level for asterisks
if p < 0.0001
    stars = '****';
elseif p < 0.001
    stars = '***';
elseif p < 0.01
    stars = '**';
elseif p < 0.05
    stars = '*';
else
    stars = 'n.s.'; % Not significant
end

% Add asterisks to indicate significance on the plot
xpos = 1.5; % Position between the two boxes
ypos = max([mean_p_val_1, mean_p_val_2]); % Adjust the y-position as needed
text(xpos, ypos, stars, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
hold off;
