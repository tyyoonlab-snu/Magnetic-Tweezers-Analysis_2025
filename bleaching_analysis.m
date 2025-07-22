function bleaching_analysis_main()
    % --- Configuration ---
    % Set the base folder containing the experimental data.
    baseFolder = '';

    if contains(baseFolder, 'RFP', 'IgnoreCase', true)
        fprintf('--- Detected RFP in folder name. Using RFP parameters. ---\n');
        fluorophore_type = 'RFP'; %mCherry 
        
        % --- Trace Extraction & Analysis Parameters ---
        roiSize = 2; % Size of the square ROI for intensity measurement.
        edgeMargin = 50; % Pixels to ignore around the image border.
        
        % --- Image Pre-processing Parameters ---
        channel_to_analyze = 'Full'; % Options: 'Top', 'Bottom', 'Full'
        frames_to_average = 3:10; 
        
        % --- Spot Finding Parameters ---
        min_peak_intensity = 200; % Adjusted for RFP
        min_peak_distance = 2; % In pixels.
       
        % --- Automatic Curation & Step-Finding Parameters ---
        median_filter_order = 11; 
        step_delta_threshold = -5; 
        delta_max_peak_width = 5;  
        recovery_delta_threshold = 7;   
        max_deviation_from_median = 25; 
        
        % --- Final Intensity Check Parameters ---
        final_intensity_threshold = 160; 
        min_trace_length_for_bleach_check = 200;
        final_bleach_check_window_fraction = 0.2; 
        
        % --- Display Parameters for Spot Detection Image ---
        display_min_intensity = 50;
        display_max_intensity = 270;
    else
        fprintf('--- Using default GFP parameters. ---\n');
        fluorophore_type = 'GFP';
        
        % --- Trace Extraction & Analysis Parameters ---
        roiSize = 2;
        edgeMargin = 50;
        
        % --- Image Pre-processing Parameters ---
        channel_to_analyze = 'Full';
        frames_to_average = 3:10; 
        
        % --- Spot Finding Parameters ---
        min_peak_intensity = 190; 
        min_peak_distance = 2; % In pixels.        
        
        % --- Automatic Curation & Step-Finding Parameters ---
        median_filter_order = 11; 
        step_delta_threshold = -5;
        delta_max_peak_width = 5;
        recovery_delta_threshold = 7;
        max_deviation_from_median = 25;     

        % --- Final Intensity Check Parameters ---
        final_intensity_threshold = 160;
        min_trace_length_for_bleach_check = 200;
        final_bleach_check_window_fraction = 0.2;
        
        % --- Display Parameters for Spot Detection Image ---
        display_min_intensity = 60; 
        display_max_intensity = 270;
    end
    
    % --- Initialization ---
    if ~isfolder(baseFolder)
        error('The specified base folder does not exist: %s', baseFolder);
    end
    % Initialize structs for storing results. Fields related to histogram analysis removed.
    allKeptResults = struct('ImageSeries', {}, 'SpotID', {}, 'ROICenter', {}, 'IntensityTrace', {}, 'BleachingSteps', {}, 'StepIndices', {}, 'DTrace', {}, 'PotentialStepIndices', {});
    allDiscardedResults = struct('ImageSeries', {}, 'SpotID', {}, 'ROICenter', {}, 'IntensityTrace', {}, 'Reason', {}, 'DiscardPoint', {});
    
    % --- Main Processing Loop ---
    fprintf('\nProcessing Base Folder: %s\n', baseFolder);
    
    spotFolders = dir(fullfile(baseFolder, '**', 'RawData_InWell*'));
    spotFolders = spotFolders([spotFolders.isdir]); 
    if isempty(spotFolders)
        fprintf('  -> No ''RawData_InWell*'' folders found within the specified directory or its subdirectories.\n');
        return;
    end
    
    for j = 1:length(spotFolders)
        seriesName = spotFolders(j).name;
        seriesFolder = spotFolders(j).folder;
        seriesPath = fullfile(seriesFolder, seriesName);
        
        relativeSeriesPath = strrep(seriesPath, baseFolder, '');
        if startsWith(relativeSeriesPath, filesep)
            relativeSeriesPath = relativeSeriesPath(2:end);
        end
        
        fprintf('  -> Analyzing Image Series: %s\n', relativeSeriesPath);
        try
            % Call to analyze_image_series updated to remove histogram and debug plot parameters.
            [keptTraces, discardedTraces] = analyze_image_series(seriesPath, channel_to_analyze, frames_to_average, ...
                min_peak_intensity, min_peak_distance, roiSize, edgeMargin, display_min_intensity, display_max_intensity, ...
                median_filter_order, step_delta_threshold, delta_max_peak_width, ...
                recovery_delta_threshold, max_deviation_from_median, final_intensity_threshold, ...
                min_trace_length_for_bleach_check, final_bleach_check_window_fraction, fluorophore_type);
            
            if ~isempty(keptTraces)
                for k = 1:length(keptTraces)
                    keptTraces(k).ImageSeries = relativeSeriesPath;
                end
                allKeptResults = [allKeptResults, keptTraces];
            end
            if ~isempty(discardedTraces)
                 for k = 1:length(discardedTraces)
                    discardedTraces(k).ImageSeries = relativeSeriesPath;
                end
                allDiscardedResults = [allDiscardedResults, discardedTraces];
            end
            
            fprintf('    - Kept %d traces, Discarded %d traces from this series.\n', length(keptTraces), length(discardedTraces));
            
        catch ME
            fprintf('    - ERROR processing series %s: %s\n', relativeSeriesPath, ME.message);
            fprintf('    - Error in function %s, line %d\n', ME.stack(1).name, ME.stack(1).line);
        end
    end
    
    % --- Results Visualization ---
    if ~isempty(allKeptResults) || ~isempty(allDiscardedResults)
        fprintf('\n--- Analysis Complete ---\n');
        % Generate the overall summary plot. All debugging plots have been removed.
        generate_summary_plot(allKeptResults);
    else
        fprintf('\nNo data was processed successfully.\n');
    end
end

% ------------------------------------------------------------------------
%                        HELPER FUNCTIONS
% ------------------------------------------------------------------------

function [keptResults, discardedResults] = analyze_image_series(seriesPath, channel, avgFrames, minIntensity, minDistance, roiSize, edgeMargin, dispMin, dispMax, filterOrder, stepThresh, deltaPeakWidth, recoveryThresh, deviationThresh, finalIntensityThresh, minLengthBleach, fractionBleach, fluorophore_type)
    
    % --- Initialize result structs for this series ---
    keptResults = struct('SpotID', {}, 'ROICenter', {}, 'IntensityTrace', {}, 'BleachingSteps', {}, 'StepIndices', {}, 'DTrace', {}, 'PotentialStepIndices', {});
    discardedResults = struct('SpotID', {}, 'ROICenter', {}, 'IntensityTrace', {}, 'Reason', {}, 'DiscardPoint', {});
    
    % --- Step 1: Load a sequence of Single-Frame TIFFs with Numerical Sorting ---
    imageFiles_tif = dir(fullfile(seriesPath, '*.tif'));
    imageFiles_tiff = dir(fullfile(seriesPath, '*.tiff'));
    imageFiles = [imageFiles_tif; imageFiles_tiff];
    if isempty(imageFiles)
        warning('No TIFF images found in: %s', seriesPath);
        return;
    end
    
    frameNumbers = zeros(numel(imageFiles), 1);
    for i = 1:numel(imageFiles)
        numeric_part = regexp(imageFiles(i).name, '\d+', 'match');
        if ~isempty(numeric_part)
            frameNumbers(i) = str2double(numeric_part{end});
        else
            frameNumbers(i) = i;
        end
    end
    
    [~, sorted_indices] = sort(frameNumbers);
    sortedImageFiles = imageFiles(sorted_indices);
    
    info = imfinfo(fullfile(seriesPath, sortedImageFiles(1).name));
    numImages = numel(sortedImageFiles);
    
    fullHeight = info(1).Height;
    imgWidth = info(1).Width;
    imgHeight = fullHeight;
    y_offset = 0;
    if strcmpi(channel, 'Top')
        imgHeight = floor(fullHeight / 2);
        y_offset = imgHeight;
    elseif strcmpi(channel, 'Bottom')
        imgHeight = floor(fullHeight / 2);
    end
    
    imgStack = zeros(imgHeight, imgWidth, numImages, 'uint16');
    for k = 1:numImages
        filePath = fullfile(seriesPath, sortedImageFiles(k).name);
        full_img = imread(filePath);
        imgStack(:,:,k) = full_img(y_offset+1 : y_offset+imgHeight, :);
    end
    
    % --- Step 2 & 3: Find and Refine Spot Locations ---
    avgImg = double(mean(imgStack(:,:,avgFrames), 3));
    peaks = imregionalmax(avgImg, 8);
    peaks(avgImg < minIntensity) = 0;
    [y_coords, x_coords] = find(peaks);
    initial_centroids = [x_coords, y_coords];
    if isempty(initial_centroids); return; end
    intensities = avgImg(sub2ind(size(avgImg), y_coords, x_coords));
    [~, sortIdx] = sort(intensities, 'descend');
    sorted_centroids = initial_centroids(sortIdx, :);
    candidate_centroids = [];
    while ~isempty(sorted_centroids)
        candidate_centroids = [candidate_centroids; sorted_centroids(1,:)];
        distances = pdist2(sorted_centroids(1,:), sorted_centroids);
        sorted_centroids(distances < minDistance, :) = [];
    end
    if edgeMargin > 0
        validIndices = candidate_centroids(:,1) > edgeMargin & candidate_centroids(:,1) < (imgWidth - edgeMargin) & ...
                       candidate_centroids(:,2) > edgeMargin & candidate_centroids(:,2) < (imgHeight - edgeMargin);
        candidate_centroids = candidate_centroids(validIndices, :);
    end
    if isempty(candidate_centroids); return; end
    
    fit_roi_half_size = floor(roiSize / 2);
    refined_centroids = zeros(size(candidate_centroids));
    spots_to_keep = false(size(candidate_centroids,1), 1);
    for s = 1:size(candidate_centroids, 1)
        x_cen_int = round(candidate_centroids(s,1));
        y_cen_int = round(candidate_centroids(s,2));
        x1 = max(1, x_cen_int - fit_roi_half_size);
        x2 = min(imgWidth, x_cen_int + fit_roi_half_size);
        y1 = max(1, y_cen_int - fit_roi_half_size);
        y2 = min(imgHeight, y_cen_int + fit_roi_half_size);
        fit_roi = avgImg(y1:y2, x1:x2);
        [fit_params, was_fit_successful] = fit_2d_gaussian(fit_roi);
        if was_fit_successful
            refined_centroids(s, 1) = x1 + fit_params(2) - 1;
            refined_centroids(s, 2) = y1 + fit_params(3) - 1;
            spots_to_keep(s) = true;
        end
    end
    centroids = refined_centroids(spots_to_keep, :);
    numSpots = size(centroids, 1);
    if numSpots == 0; return; end
    
    % --- Step 4: Extract Intensity Traces ---
    intensityTraces = zeros(numSpots, numImages);
    halfRoi = floor(roiSize / 2);
    for s = 1:numSpots
        x_cen = centroids(s,1);
        y_cen = centroids(s,2);
        x1 = max(1, round(x_cen - halfRoi));
        x2 = min(imgWidth, round(x_cen + halfRoi));
        y1 = max(1, round(y_cen - halfRoi));
        y2 = min(imgHeight, round(y_cen + halfRoi));
        for k = 1:numImages
            signal_roi = imgStack(y1:y2, x1:x2, k);
            intensityTraces(s, k) = mean(signal_roi(:));
        end
    end
    
    % --- Step 5: Automatic Curation and Analysis ---
    for s = 1:numSpots
        trace = intensityTraces(s, :);
        
        % --- Apply Curation Filters ---
        [is_recovered, reason_mono, recovery_idx] = check_for_recovery_by_delta(trace, recoveryThresh, filterOrder);
        [is_noisy, reason_noise, noise_idx] = check_for_noise(trace, deviationThresh, filterOrder);
        [is_bleached, reason_bleach] = check_final_intensity(trace, finalIntensityThresh, minLengthBleach, fractionBleach);
        
        if ~is_recovered && ~is_noisy && is_bleached
            newEntry = struct(); % Initialize struct to prevent field mismatch
            
            % Step finding now only uses the delta method.
            [numSteps, step_indices, delta_data] = find_steps_by_delta(trace, stepThresh, filterOrder, fluorophore_type, deltaPeakWidth);
            
            newEntry.SpotID = s;
            newEntry.ROICenter = centroids(s, :);
            newEntry.IntensityTrace = trace;
            newEntry.BleachingSteps = numSteps;
            newEntry.StepIndices = step_indices;
            
            if ~isempty(delta_data)
                newEntry.DTrace = delta_data.d_trace;
                newEntry.PotentialStepIndices = delta_data.potential_step_indices;
            else
                newEntry.DTrace = [];
                newEntry.PotentialStepIndices = [];
            end
            
            keptResults(end+1) = newEntry;
        else
            newDiscard = struct(); % Initialize struct to prevent field mismatch
            newDiscard.SpotID = s;
            newDiscard.ROICenter = centroids(s, :);
            newDiscard.IntensityTrace = trace;
            if ~is_bleached
                newDiscard.Reason = reason_bleach;
                newDiscard.DiscardPoint = [];
            elseif is_recovered
                newDiscard.Reason = reason_mono;
                newDiscard.DiscardPoint = recovery_idx;
            else
                newDiscard.Reason = reason_noise;
                newDiscard.DiscardPoint = noise_idx;
            end
            discardedResults(end+1) = newDiscard;
        end
    end
    
    % --- Step 6: Save Summary of Step Counts ---
    if ~isempty(keptResults)
        wellFolderPath = fileparts(seriesPath);
        
        % Simplified folder naming, as only one method is used.
        summaryFolderName = sprintf('analysis_summary-intensity_%d-step_delta_%d-width_%d', minIntensity, abs(stepThresh), deltaPeakWidth);
        
        summaryFolder = fullfile(wellFolderPath, summaryFolderName);
        if ~exist(summaryFolder, 'dir')
           mkdir(summaryFolder)
        end
        
        step_counts = [keptResults.BleachingSteps];
        max_step = max([step_counts, 0]);
        step_bins = 0:max_step;
        counts_per_step = histcounts(step_counts, [step_bins, max_step + 1] - 0.5);
        
        output_table = table(step_bins', counts_per_step', 'VariableNames', {'NumberOfSteps', 'TraceCount'});
        [~, seriesNameOnly, ~] = fileparts(seriesPath);
        csvFileName = [seriesNameOnly, '.csv'];
        summary_filename = fullfile(summaryFolder, csvFileName);
        writetable(output_table, summary_filename);
    end
end

function [isRecovered, reason, discard_idx] = check_for_recovery_by_delta(trace, recovery_delta, filter_order)
    isRecovered = false;
    reason = 'Good';
    discard_idx = [];
    if length(trace) < 10; return; end
    
    filtered_trace = medfilt1(trace, filter_order, 'truncate');
    deltas = diff(filtered_trace);
    
    recovery_point = find(deltas > recovery_delta, 1, 'first');
    
    if ~isempty(recovery_point)
        isRecovered = true;
        reason = 'Recovery';
        discard_idx = recovery_point + 1;
    end
end

function [isNoisy, reason, discard_idx] = check_for_noise(trace, deviation_threshold, filter_order)
    isNoisy = false;
    reason = 'Good';
    discard_idx = [];
    if length(trace) < 10; return; end
    
    if mod(filter_order, 2) == 0; filter_order = filter_order + 1; end
    
    filtered_trace = medfilt1(trace, filter_order, 'truncate');
    deviations = abs(trace - filtered_trace);
    noise_point = find(deviations > deviation_threshold, 1, 'first');
    
    if ~isempty(noise_point)
        isNoisy = true;
        reason = 'High Noise';
        discard_idx = noise_point;
    end
end

function [isBleached, reason] = check_final_intensity(trace, threshold, min_length, window_fraction)
    isBleached = true;
    reason = 'Good';
    trace_length = length(trace);
    
    if trace_length < min_length
        isBleached = false;
        reason = 'Too Short';
        return;
    end
    
    start_frame = round(trace_length * (1 - window_fraction)) + 1;
    if start_frame > trace_length
        start_frame = trace_length;
    end
    
    final_avg = mean(trace(start_frame:end));
    
    if final_avg >= threshold
        isBleached = false;
        reason = 'Not Bleached';
    end
end

function [numSteps, step_indices, delta_data] = find_steps_by_delta(trace, step_delta, filter_order, fluorophore_type, max_peak_width)
    
    % --- Initializations ---
    delta_data = struct('d_trace', [], 'potential_step_indices', []);
    numSteps = 0;
    step_indices = [];
    % --- Configuration for Step Finding ---
    min_level_intensity_diff = 10; % The mean intensity of consecutive levels must differ by at least this value.
    min_frames_between_steps = 10; % Minimum distance between steps in frames.
    
    % --- Initial Checks ---
    if isempty(trace) || length(trace) < min_frames_between_steps
        return;
    end
    
    % --- Filter the trace to reduce noise ---
    if mod(filter_order, 2) == 0; filter_order = filter_order + 1; end
    filtered_trace = medfilt1(trace, filter_order, 'truncate');
    
    % --- Find sharp drops using the first derivative ---
    d_trace = -diff(filtered_trace);
    delta_data.d_trace = d_trace;
    
    min_drop_magnitude = abs(step_delta);
    
    % Find locations of sharp drops meeting magnitude AND width requirements.
    [~, potential_step_indices] = findpeaks(d_trace, ...
        'MinPeakHeight', min_drop_magnitude, ...
        'MaxPeakWidth', max_peak_width);
    
    potential_step_indices = potential_step_indices + 1;
    delta_data.potential_step_indices = potential_step_indices;
    if isempty(potential_step_indices)
        return;
    end
    
    % --- Filter the potential steps to find the true, distinct steps ---
    accepted_steps = potential_step_indices(1);
    
    for i = 2:length(potential_step_indices)
        current_potential_step = potential_step_indices(i);
        last_accepted_step = accepted_steps(end);
        
        if (current_potential_step - last_accepted_step) < min_frames_between_steps
            continue;
        end
        
        window_size = 5;
        last_level_end_idx = min(length(filtered_trace), last_accepted_step + window_size);
        mean_intensity_last_level = mean(filtered_trace(last_accepted_step:last_level_end_idx));
        
        current_level_end_idx = min(length(filtered_trace), current_potential_step + window_size);
        mean_intensity_current_level = mean(filtered_trace(current_potential_step:current_level_end_idx));
        
        if (mean_intensity_last_level - mean_intensity_current_level) >= min_level_intensity_diff
            accepted_steps = [accepted_steps, current_potential_step];
        end
    end
    
    step_indices = accepted_steps;
    
    if strcmpi(fluorophore_type, 'RFP')
        max_step=4;
    else 
        max_step=4;
    end
    
    if length(step_indices) > max_step
        accepted_step_starts = step_indices - 1;
        accepted_step_starts(accepted_step_starts < 1) = 1;
        
        valid_delta_indices = accepted_step_starts(accepted_step_starts <= length(d_trace));
        accepted_delta_values = d_trace(valid_delta_indices);
        
        [~, sorted_idx] = sort(accepted_delta_values, 'descend');
        
        top_indices = step_indices(sorted_idx(1:max_step));
        step_indices = sort(top_indices);
    end
    
    numSteps = length(step_indices);
end

function [fit_params, success] = fit_2d_gaussian(roi)
    [height, width] = size(roi);
    [X, Y] = meshgrid(1:width, 1:height);
    gaussian_model = @(p, xy) p(1) * exp(-(((xy(:,1)-p(2)).^2)/(2*p(4)^2) + ((xy(:,2)-p(3)).^2)/(2*p(4)^2))) + p(5);
    [max_val, max_idx] = max(roi(:));
    [y0_guess, x0_guess] = ind2sub(size(roi), max_idx);
    bg_guess = min(roi(:));
    amp_guess = max_val - bg_guess;
    sigma_guess = 1.5;
    p0 = [amp_guess, x0_guess, y0_guess, sigma_guess, bg_guess];
    lb = [0, 1, 1, 0.1, 0];
    ub = [double(intmax('uint16')), width, height, 5, double(intmax('uint16'))];
    xy_data = [X(:), Y(:)];
    z_data = roi(:);
    options = optimoptions('lsqcurvefit', 'Display', 'none');
    try
        fit_params = lsqcurvefit(gaussian_model, p0, xy_data, z_data, lb, ub, options);
        success = true;
    catch
        fit_params = p0;
        success = false;
    end
end

function generate_summary_plot(allKeptResults)
    if isempty(allKeptResults)
        fprintf('No kept traces to summarize.\n');
        return;
    end
    conditionData = containers.Map('KeyType', 'char', 'ValueType', 'any');
    maxStepsFound = 0;
    
    % Aggregate data by condition and replicate
    for i = 1:length(allKeptResults)
        result = allKeptResults(i);
        pathParts = strsplit(result.ImageSeries, filesep);
        conditionName = pathParts{1};
        replicateName = pathParts{2};
        
        if ~isKey(conditionData, conditionName)
            conditionData(conditionName) = containers.Map('KeyType', 'char', 'ValueType', 'any');
        end
        replicateMap = conditionData(conditionName);
        
        if ~isKey(replicateMap, replicateName)
            % Initialize with a buffer size, will be trimmed later
            replicateMap(replicateName) = zeros(1, 10);
        end
        
        counts = replicateMap(replicateName);
        step = result.BleachingSteps;
        
        % Ensure the counts array is large enough
        if step + 1 > length(counts)
            counts(step + 1) = 0;
        end
        
        if step >= 0 % Include 0-step events
            counts(step + 1) = counts(step + 1) + 1;
            replicateMap(replicateName) = counts;
            conditionData(conditionName) = replicateMap;
            maxStepsFound = max(maxStepsFound, step);
        end
    end
    
    if isempty(keys(conditionData))
        fprintf('No steps were found in any kept traces.\n');
        return;
    end
    
    conditionNames = keys(conditionData);
    numConditions = length(conditionNames);
    
    meanData = zeros(numConditions, maxStepsFound + 1);
    stdData = zeros(numConditions, maxStepsFound + 1);
    legendLabels = cell(1, numConditions);
    
    for i = 1:numConditions
        conditionName = conditionNames{i};
        replicateMap = conditionData(conditionName);
        replicateNames = keys(replicateMap);
        
        % Pre-allocate with the correct size
        replicateMatrix = zeros(length(replicateNames), maxStepsFound + 1);
        
        for j = 1:length(replicateNames)
            counts = replicateMap(replicateNames{j});
            len = length(counts);
            replicateMatrix(j, 1:min(len, maxStepsFound + 1)) = counts(1:min(len, maxStepsFound + 1));
        end
        
        totalSpots = sum(replicateMatrix(:));
        legendLabels{i} = sprintf('%s (n=%d)', strrep(conditionName, '_', '\_'), totalSpots);
        
        meanData(i, :) = mean(replicateMatrix, 1);
        stdData(i, :) = std(replicateMatrix, 0, 1);
    end
    
    % --- Plot Observed Step Distribution ---
    figure('Name', 'Summary of Observed Bleaching Steps', 'NumberTitle', 'off', 'Position', [200, 200, 900, 600]);
    
    % We plot from step 1 onwards, so we index from column 2
    if maxStepsFound > 0
        b = bar(meanData(:, 2:end)', 'grouped');
        hold on;
        
        numBarsPerGroup = size(meanData, 1);
        numStepCategories = size(meanData, 2) - 1;
        
        for i = 1:numBarsPerGroup
            x_coords = b(i).XEndPoints;
            y_data = b(i).YData;
            std_data_for_series = stdData(i, 2:end);
            errorbar(x_coords, y_data, std_data_for_series, 'k', 'linestyle', 'none');
        end
        hold off;
        
        ax = gca;
        ax.XTickLabel = 1:maxStepsFound;
        xlabel('Number of Bleaching Steps');
    else
        bar([]); % Create an empty plot if only 0-step traces were found
        xlabel('Number of Bleaching Steps (None > 0 Found)');
    end
    
    title('Distribution of Observed Bleaching Steps by Condition');
    ylabel('Average Count per Replicate');
    legend(legendLabels, 'Location', 'northeast');
    grid on;
    fprintf('\nSummary plot of observed steps has been generated.\n');
end
