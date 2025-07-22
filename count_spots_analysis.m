
function spot_counting_main()
    % --- Configuration ---
    % Base folder containing the well subfolders (A3, A4, etc.).
    baseFolder = '';
    
    % --- Debugging ---
    % Set to 1 to display an image for EACH analyzed image showing detected spots.
    % Set to 0 to run without generating plots (except for one example).
    showDebugPlot = 1; % Set to 1 to generate an example detection figure.
    
    % --- Spot Finding Parameters ---
    % The minimum intensity a pixel must have to be considered a peak.
    min_peak_intensity = 130; 
    % The minimum distance (in pixels) between two separate spots. 
    min_peak_distance = 2; 
    % Number of pixels to ignore around the image border to avoid edge artifacts.
    edgeMargin = 50; 
    
    % --- Parameters for Lab-Specific Peak Finder ---
    % Half-size of the search area for finding local maxima (e.g., 4 = 9x9 box).
    localarea_size = 4;
    % Set to true to apply a Gaussian filter before peak finding.
    use_gauss_filter = true;
    % Sigma value for the Gaussian filter if used.
    gauss_filter_sigma = 1.5;
    
    % --- Display Parameters for Spot Detection Image ---
    % These control the contrast of the debug plot image.
    display_min_intensity = 50; 
    display_max_intensity = 270;
    
    % --- Initialization ---
    if ~isfolder(baseFolder)
        error('The specified base folder does not exist: %s', baseFolder);
    end
    
    % Structure to store the final statistics.
    allResults = struct('Well', {}, 'AverageSpotCount', {}, 'StdDevSpotCount', {}, 'AverageSummedIntensity', {}, 'StdDevSummedIntensity', {});
    
    % --- Main Processing Loop ---
    fprintf('\nProcessing Base Folder: %s\n', baseFolder);
    
    % Find all subdirectories in the base folder.
    subFolders = dir(baseFolder);
    allSubFolders = subFolders([subFolders.isdir]);
    % Filter out '.' and '..' directories
    allSubFolders = allSubFolders(~ismember({allSubFolders.name}, {'.', '..'}));

    % Further filter for folders with names like A3, C5, etc. (letter followed by number)
    folderNames = {allSubFolders.name};
    % Regular expression to match a name starting with a letter, followed by one or more digits.
    % e.g., 'A3', 'C10', 'D8' will match. 'A', 'Data', 'A3a' will not.
    matchingFoldersIdx = ~cellfun('isempty', regexp(folderNames, '^[A-Za-z]\d+$'));
    wellFolders = allSubFolders(matchingFoldersIdx);
    
    if isempty(wellFolders)
        fprintf('  -> No well folders matching the pattern (e.g., A3, C5) found in the base directory.\n');
        return;
    end
    
    examplePlotShown = false;
    for j = 1:length(wellFolders)
        wellName = wellFolders(j).name;
        wellPath = fullfile(baseFolder, wellName);
        
        fprintf('  -> Analyzing Well: %s\n', wellName);
        
        try
            % Determine if a plot should be shown for the current well.
            showPlotForThisWell = showDebugPlot && ~examplePlotShown;
            
            % This function analyzes all images in a well and returns stats for counts and summed intensity.
            [avgCount, stdDevCount, avgSummedIntensity, stdDevSummedIntensity, plotWasShown] = analyze_well_folder(wellPath, min_peak_intensity, min_peak_distance, edgeMargin, showPlotForThisWell, display_min_intensity, display_max_intensity, localarea_size, use_gauss_filter, gauss_filter_sigma);
            
            if plotWasShown
                examplePlotShown = true;
            end

            if ~isnan(avgCount)
                % Store the result
                newEntry.Well = wellName;
                newEntry.AverageSpotCount = avgCount;
                newEntry.StdDevSpotCount = stdDevCount;
                newEntry.AverageSummedIntensity = avgSummedIntensity;
                newEntry.StdDevSummedIntensity = stdDevSummedIntensity;
                allResults(end+1) = newEntry;

                fprintf('    - Avg Count: %.2f (Std: %.2f) | Avg Summed Intensity: %.2e (Std: %.2e)\n', avgCount, stdDevCount, avgSummedIntensity, stdDevSummedIntensity);
            else
                fprintf('    - No valid images found or processed for this well.\n');
            end
            
        catch ME
            fprintf('    - ERROR processing well %s: %s\n', wellName, ME.message);
            fprintf('    - Error in function %s, line %d\n', ME.stack(1).name, ME.stack(1).line);
        end
    end
    
    % --- Save Results and Generate Plot ---
    if ~isempty(allResults)
        fprintf('\n--- Analysis Complete ---\n');
        
        % Convert the results structure to a table for easy saving.
        resultsTable = struct2table(allResults);
        
        % Sort the table by Well name
        resultsTable = sortrows(resultsTable, {'Well'});
        
        % Save the summary table to a CSV file in the base folder.
        outputFileName = fullfile(baseFolder, 'spot_statistics.csv');
        writetable(resultsTable, outputFileName);
        fprintf('Summary of spot statistics saved to: %s\n', outputFileName);
        
        % Generate and display the final summary plots
        generate_count_summary_plot(resultsTable);
        generate_summed_intensity_summary_plot(resultsTable);
        fprintf('Summary plots have been generated.\n');

    else
        fprintf('\nNo data was processed successfully.\n');
    end
end

function [avgCount, stdDevCount, avgSummedIntensity, stdDevSummedIntensity, plotWasShown] = analyze_well_folder(wellPath, minIntensity, minDistance, edgeMargin, showPlot, dispMin, dispMax, localarea_size, use_gauss_filter, gauss_filter_sigma)
    % This function finds all pre-averaged images in a well folder, counts
    % the spots, sums their intensities, and returns statistics.
    
    avgCount = NaN;
    stdDevCount = NaN;
    avgSummedIntensity = NaN;
    stdDevSummedIntensity = NaN;
    plotWasShown = false;
    
    % Find all files ending with '_AVG.tif' or '_AVG.tiff'
    imageFiles_tif = dir(fullfile(wellPath, '*_AVG.tif'));
    imageFiles_tiff = dir(fullfile(wellPath, '*_AVG.tiff'));
    imageFiles = [imageFiles_tif; imageFiles_tiff];
    
    if isempty(imageFiles)
        warning('No ''*_AVG.tif'' images found in: %s', wellPath);
        return;
    end
    
    [imgHeight, imgWidth] = size(imread(fullfile(wellPath, imageFiles(1).name)));
    all_spot_counts = [];
    all_summed_intensities_per_image = []; % To store the summed intensity for each image
    
    % Process each image in the well folder
    for i = 1:length(imageFiles)
        filePath = fullfile(wellPath, imageFiles(i).name);
        avgImg = double(imread(filePath));
        
        % Find peaks using the custom lab-specific method
        [x_coords, y_coords] = find_peaks_lab_method(avgImg, localarea_size, minIntensity, use_gauss_filter, gauss_filter_sigma);
        
        initial_centroids = [x_coords, y_coords];
        if isempty(initial_centroids)
            all_spot_counts(end+1) = 0;
            all_summed_intensities_per_image(end+1) = 0; % Record zero intensity if no spots
            continue; % No spots found in this image
        end
        
        % --- Refine Spot Locations ---
        
        % Filter out peaks that are too close to each other, keeping the brightest
        intensities = avgImg(sub2ind(size(avgImg), y_coords, x_coords));
        [~, sortIdx] = sort(intensities, 'descend');
        sorted_centroids = initial_centroids(sortIdx, :);
        
        candidate_centroids = [];
        while ~isempty(sorted_centroids)
            candidate_centroids = [candidate_centroids; sorted_centroids(1,:)];
            distances = pdist2(sorted_centroids(1,:), sorted_centroids);
            sorted_centroids(distances < minDistance, :) = [];
        end
        
        % Remove spots too close to the edge
        if edgeMargin > 0
            validIndices = candidate_centroids(:,1) > edgeMargin & candidate_centroids(:,1) < (imgWidth - edgeMargin) & ...
                           candidate_centroids(:,2) > edgeMargin & candidate_centroids(:,2) < (imgHeight - edgeMargin);
            final_centroids = candidate_centroids(validIndices, :);
        else
            final_centroids = candidate_centroids;
        end
        
        spotCount = size(final_centroids, 1);
        all_spot_counts(end+1) = spotCount;
        
        % --- Measure and Store SUM of Intensities ---
        if spotCount > 0
            % Get the linear indices of the final centroids
            linear_indices = sub2ind(size(avgImg), round(final_centroids(:,2)), round(final_centroids(:,1)));
            % Extract the intensity value at each spot's location
            current_intensities = avgImg(linear_indices);
            % Sum the intensities for this image and add to our master list for the well
            all_summed_intensities_per_image(end+1) = sum(current_intensities);
        else
            all_summed_intensities_per_image(end+1) = 0; % Record zero intensity if no spots
        end
        
        % --- Display Debug Plot (Optional) ---
        if showPlot && ~plotWasShown && spotCount > 0
            figure('Name', ['Example Spot Detection: ' strrep(imageFiles(i).name, '_', '\_')], 'NumberTitle', 'off');
            imshow(avgImg, [dispMin, dispMax]); 
            hold on;
            viscircles(final_centroids, repmat(minDistance, spotCount, 1), 'Color', 'g', 'LineWidth', 0.5);
            hold off;
            title(sprintf('%d spots found', spotCount));
            plotWasShown = true; % Mark that the plot has been shown
        end
    end
    
    % Calculate statistics for the well
    if ~isempty(all_spot_counts)
        avgCount = mean(all_spot_counts);
        stdDevCount = std(all_spot_counts);
    end
    if ~isempty(all_summed_intensities_per_image)
        avgSummedIntensity = mean(all_summed_intensities_per_image);
        stdDevSummedIntensity = std(all_summed_intensities_per_image);
    end
end

function [x_coords, y_coords] = find_peaks_lab_method(frame, localarea_size, min_intensity, use_gauss_filter, sigma)
    % Finds local maxima using the lab's original iterative method.
    % 1. Optionally applies a Gaussian filter.
    % 2. Iterates through each pixel.
    % 3. Checks if the pixel is the maximum in its local neighborhood.
    % 4. Checks if the pixel intensity is above a minimum threshold.

    % Apply Gaussian filter if requested
    if use_gauss_filter
        frame = imgaussfilt(frame, sigma);
    end
    
    [Xsize, Ysize] = size(frame);
    border = localarea_size;
    
    x_coords = [];
    y_coords = [];
    
    % Iterate through the image, avoiding the borders
    for i = (border + 1):(Xsize - border)
        for j = (border + 1):(Ysize - border)
            
            current_pixel_intensity = frame(i, j);
            
            % Condition 1: Check if intensity is above the minimum threshold
            if current_pixel_intensity > min_intensity
                
                % Define the local area to check for a maximum
                subframe = frame(i - localarea_size : i + localarea_size, ...
                                 j - localarea_size : j + localarea_size);
                
                % Condition 2: Check if the current pixel is the local maximum
                if current_pixel_intensity >= max(subframe(:))
                    % This is a peak. Add its coordinates.
                    % Note: The check is >= to handle plateaus, where one of
                    % the plateau pixels will be chosen.
                    x_coords(end+1, 1) = j;
                    y_coords(end+1, 1) = i;
                end
            end
        end
    end
end

function generate_count_summary_plot(resultsTable)
    % This function creates a bar chart of the average spot counts per well
    % with standard deviation error bars.
    
    if isempty(resultsTable) || ~ismember('AverageSpotCount', resultsTable.Properties.VariableNames)
        disp('No count data to plot.');
        return;
    end
    
    figure('Name', 'Summary of Spot Counts per Well', 'NumberTitle', 'off', 'Position', [200, 200, 800, 600]);
    
    % Extract data for plotting
    wellNames = resultsTable.Well;
    avgCounts = resultsTable.AverageSpotCount;
    stdDevs = resultsTable.StdDevSpotCount;
    
    % Create the bar chart
    b = bar(avgCounts);
    
    % Add error bars
    hold on;
    x_coords = b.XEndPoints;
    errorbar(x_coords, avgCounts, stdDevs, 'k', 'linestyle', 'none', 'LineWidth', 1.5);
    hold off;
    
    % Customize the plot
    ax = gca;
    ax.XTickLabel = wellNames;
    ax.XTickLabelRotation = 45; 
    
    title('Average Spot Count per Well');
    xlabel('Well');
    ylabel('Average Number of Spots');
    grid on;
end

function generate_summed_intensity_summary_plot(resultsTable)
    % This function creates a bar chart of the average of the summed spot
    % intensities per well with standard deviation error bars.
    
    if isempty(resultsTable) || ~ismember('AverageSummedIntensity', resultsTable.Properties.VariableNames)
        disp('No intensity data to plot.');
        return;
    end
    
    figure('Name', 'Summary of Summed Spot Intensity per Well', 'NumberTitle', 'off', 'Position', [1050, 200, 800, 600]);
    
    % Extract data for plotting
    wellNames = resultsTable.Well;
    avgIntensities = resultsTable.AverageSummedIntensity;
    stdDevs = resultsTable.StdDevSummedIntensity;
    
    % Create the bar chart
    b = bar(avgIntensities);
    
    % Add error bars
    hold on;
    x_coords = b.XEndPoints;
    errorbar(x_coords, avgIntensities, stdDevs, 'k', 'linestyle', 'none', 'LineWidth', 1.5);
    hold off;
    
    % Customize the plot
    ax = gca;
    ax.XTickLabel = wellNames;
    ax.XTickLabelRotation = 45;
    
    title('Average of Summed Spot Intensities per Well');
    xlabel('Well');
    ylabel('Average Summed Intensity (a.u.)');
    grid on;
end
