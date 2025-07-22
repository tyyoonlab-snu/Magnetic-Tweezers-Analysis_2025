%% Load files
close all;
clear; clc;

% Prompt the user to select a folder
folder = uigetdir();
if folder == 0
    error('No folder selected. Exiting script.');
end
cd(folder);

% Get Room number
RM_number = [];
while isempty(RM_number)
    RM_number = input('Room number? (1 or 2): ');
    if ~ismember(RM_number, [1, 2])
        RM_number = [];
        disp('Invalid input. Please enter 1 or 2.');
    end
end

% Define parameters based on room number
switch RM_number
    case 1
        F0 = -0.0667; A1 = 0.0113; d1 = -3.3477;
        A2 = 2.4916e-07; d2 = -1.3388;
    case 2
        F0 = 0.007229; A1 = 0.02587; d1 = -3.732736;
        A2 = 4.959e-10; d2 = -0.989119;
end

% Get recording frequency
freq = [];
while isempty(freq)
    freq = input('Recording frequency? (100 for STD / 1200 for HS): ');
    if ~ismember(freq, [100, 1200])
        freq = [];
        disp('Invalid input. Please enter 100 or 1200.');
    end
end

% Set bead number
which_bead = '';
while isempty(which_bead)
    which_bead = input('Which bead do you want to analyze: ', 's');
end

% Experimental parameters
Rbead = 1400; % Bead radius (nm)
T = 298; % Temperature (K)
correction_factor = 0.878;
pixel_size = 80; % Pixel size in nm

% Locate bead folder
bead_to_analyze = ['B', which_bead, '*'];
temp = dir(bead_to_analyze);

if isempty(temp)
    error('No matching bead folder found. Check your input.');
end

%% Initialize Variables
pths = arrayfun(@(i) {[temp(i).name, '/']}, 1:numel(temp));  % Create paths for each bead
npth = numel(pths);  % Number of beads to analyze

% Create empty cell arrays to store data
[finfo, fname, nbead, nframe, fps, Roff, ori, f, t, t2, M, R, P, F, rz, rx, ry, x, y, z, dx, dy, dz] = deal(cell(npth, 1));

% Folder paths and bead identifiers
paths = temp.folder;  % Folder name
n_bead = temp.name;  % Bead identifiers
pos_BS = find(paths == '\');  % Position of backslashes in folder paths (for date parsing)
nfile = zeros(npth, 1);  % Initialize file count for each bead

% Process each bead folder
for p = 1:npth
    % List files starting with 'r*.xls' in the current bead folder
    finfo{p} = dir([pths{p}, 'r*.xls']);
    nfile(p) = numel(finfo{p});  % Number of files starting with 'r*'
    
    % Initialize data cells for each bead
    [fname{p}, Roff{p}, ori{p}, f{p}, M{p}, R{p}, P{p}, F{p}, x{p}, y{p}, z{p}, dx{p}, dy{p}, dz{p}] = deal(cell(nfile(p), 1));
    [nbead{p}, fps{p}, nframe{p}] = deal(zeros(nfile(p), 1));
    
    % Process each file for the current bead
    for n = 1:nfile(p)
        % Display progress
        disp([int2str(n / nfile(p) * 100), '% of ', pths{p}(1:end-1), '...']);
        
        % Load the file information
        fname{p}{n} = finfo{p}(n).name;
        fname_motor = ['s', fname{p}{n}(2:end)];  % Motor file name
        
        % Check if motor data file exists
        if exist([pths{p}, fname_motor], 'file')
            % Load the motor data
            dat = dlmread([pths{p}, fname_motor]);
            t2{p}{n} = dat(:, 1);
            M{p}{n} = dat(:, 2);
            F{p}{n} = F0 + A1 * exp((-M{p}{n}) / d1) + A2 * exp((-M{p}{n}) / d2);  % Force calculation
            R{p}{n} = (dat(:, 3) - floor(dat(:, 3))) * 360;  % Rotation angle
            P{p}{n} = dat(:, 4);  % Piezo position
            
            % Load bead data
            dat = dlmread([pths{p}, fname{p}{n}]);
            nframe{p}(n) = size(dat, 1);  % Number of frames
            tmp = dlmread([pths{p}, 'c', fname{p}{n}(2:4), '.fps']);
            fps{p}(n) = tmp(1);  % Frame rate
            Roff{p}{n} = tmp(2, :);
            ori{p}{n} = tmp(3, :);
            f{p}{n} = dat(:, 1);
            dat = dat(:, 2:end);
            t{p}{n} = f{p}{n} / fps{p}(n);  % Time vector
            nbead{p}(n) = size(dat, 2) / 3 - 1;  % Number of beads
            
            % Subtract XY offset
            dat(:, [1:3:end, 2:3:end]) = dat(:, [1:3:end, 2:3:end]) - repmat(mean(dat(31:60, [1:3:end, 2:3:end]), 1), [nframe{p}(n), 1]);
            
            % Store bead coordinates and calculate displacement
            rx{p}{n} = dat(:, 1) * pixel_size;
            ry{p}{n} = dat(:, 2) * pixel_size;
            rz{p}{n} = dat(:, 3);
            x{p}{n} = dat(:, 4:3:end) * pixel_size;
            y{p}{n} = dat(:, 5:3:end) * pixel_size;
            z{p}{n} = dat(:, 6:3:end);  % Store all traces
            dx{p}{n} = (x{p}{n} - repmat(rx{p}{n}, [1, nbead{p}(n)]));
            dy{p}{n} = (y{p}{n} - repmat(ry{p}{n}, [1, nbead{p}(n)]));
            dz{p}{n} = (z{p}{n} - repmat(rz{p}{n}, [1, nbead{p}(n)])) * correction_factor;
            
            % Synchronize motor data with bead data
            M{p}{n} = interp1(t2{p}{n}, M{p}{n}, t{p}{n});
            F{p}{n} = interp1(t2{p}{n}, F{p}{n}, t{p}{n});
            R{p}{n} = interp1(t2{p}{n}, R{p}{n}, t{p}{n});
            P{p}{n} = interp1(t2{p}{n}, P{p}{n}, t{p}{n});
        else
            % If motor data file does not exist, skip the current file
            warning('File %s does not exist. Skipping...', [pths{p}, fname_motor]);
            nfile(p) = nfile(p) - 1;  % Decrement file count
        end
    end
end

%% Combine Data
Fdat = [];
dxdat = [];
dydat = [];
dzdat = [];
rzdat = [];
Mdat = [];

% Concatenate data from all beads and files
for p = 1:npth
    for n = 1:nfile(p)
        Fdat = [Fdat; F{p}{n}];
        dxdat = [dxdat; dx{p}{n}];
        dydat = [dydat; dy{p}{n}];
        dzdat = [dzdat; dz{p}{n}];
        rzdat = [rzdat; rz{p}{n}];
        Mdat = [Mdat; M{p}{n}];
    end
end

r_Fdat = round(Fdat, 1);  % Round the force data to 1 decimal place


%% Data Processing
number = 1:numel(Fdat);  % Time index for plotting
number = number / freq;  % Convert to time in seconds

% Apply median filter to dzdat for smoothing
dz_med = medfilt1(dzdat, freq / 10);

%% Plotting
figure(1);

% Plot dzdat and dz_med
ax_dz = subplot(4, 1, [1, 2, 3]);
p_dz = plot(number, dzdat, 'color', [0.7 0.7 0.7]);  % Raw data in light gray
hold on;
p_dz_med = plot(number, dz_med, 'k', 'LineWidth', 0.5);  % Smoothed data in black

% Set y-axis limits and labels for dzdat
ylim([-300 520]);
ylabel('Extension (nm)');
set(gca, 'TickDir', 'in', 'FontSize', 15, 'LineWidth', 2, 'FontWeight', 'bold', 'TickLength', [0.005 0.005]);

% Force plot
ax_F = subplot(4, 1, 4);
p_F = plot(number, Fdat, 'k', 'LineWidth', 0.75);  % Force data in black

% Set force plot properties
yticks(0:10:50);
grid on;
xlabel('Time (s)');
ylabel('Force (pN)');
set(gca, 'TickDir', 'in', 'FontSize', 15, 'LineWidth', 1, 'FontWeight', 'bold', 'TickLength', [0.005 0.005]);

% Link x-axes of dzdat and force plots
linkaxes([ax_dz, ax_F], 'x');

% Set figure size and position
set(gcf, 'Position', [100, 100, 1800, 900]);

