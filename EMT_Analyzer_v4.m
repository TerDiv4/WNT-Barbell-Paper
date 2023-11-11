%% Preamble
clearvars;
close all;


% Set max number of replicates default 10

num_replicates_all = 10;

% Get the number of time points from the user
prompt = {'Enter the number of time points:'};
dlgtitle = 'Time points';
dims = [1 35];
definput = {'1'};
answerTP = inputdlg(prompt,dlgtitle,dims,definput);
num_time_points = str2double(answerTP{1});



% Prompt the user to input the number and names of treatment groups
prompt = {'Enter the number of treatment groups:'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'1'};
answerTG = inputdlg(prompt, dlgtitle, dims, definput);
num_groups = str2double(answerTG{1});
group_names = cell(num_groups,1);

% Loop over each time point and ask the user to input the number of replicates

% for j = 1:num_time_points
%     num_replicates = inputdlg(sprintf('Enter the number of replicates for time point %d', j), sprintf('Replicates for time point %d', j), 1);
%     num_replicates_all(j) = str2double(num_replicates{1});
% end

numberFiles = (num_time_points+num_replicates_all)*num_groups;

% Initialize vectors to store the data for each group
avg_dist_all = zeros(numberFiles, 1);
time_points_all = zeros(numberFiles, 1);
group_all = cell(numberFiles, 1);

%% Main Script


for i = 1:num_groups    % Loop over each treatment group

    prompt = sprintf('Enter the name of treatment group %d:', i);
    group_names{i} = char(inputdlg(prompt, dlgtitle, dims));
     
    for j = 1:num_time_points   % Loop over each time point and replicate

        % Prompt the user to input the number of replicates for the current time point
        num_replicates_current = inputdlg(sprintf('Enter the number of replicates for time point %d:', j));
        % Convert the input to a numeric value
         num_replicates_current = str2double(num_replicates_current{1});
         num_replicates_all(i,j) = num_replicates_current;

        for k = 1:num_replicates_current
           
        
            % Prompt the user to select the input image
            [filename_single, pathname_single] = uigetfile({'*.jpg;*.png;*.tif', 'Image files (*.jpg, *.png, *.tif)'}, ...
                sprintf('Select input image for time point %d, replicate %d, treatment group %s', j, k, group_names{i}), 'MultiSelect', 'off');
            cd(pathname_single);
            % Check if the user canceled the file selection dialog
            if isequal(filename_single,0) || isequal(pathname_single,0)
                disp('User canceled the file selection dialog');
                return;
            end

            % Load the selected input image
            img = imread(fullfile(pathname_single, filename_single));

            % Display the current image
            imshow(img);

            % Allow the user to select the center point of the image
            [x, y] = ginput(1);

            % Allow the user to draw a polygon around the outline
            h = drawpolygon;

            % Get the x and y coordinates of the vertices of the polygon
            xy = h.Position;

            % Compute the distance from each point on the outline to the center point
            distances = sqrt((xy(:,1)-x).^2 + (xy(:,2)-y).^2);

            % Compute the average distance
            avg_dist = mean(distances);

            % Store the average distance, time point, and treatment group for the current file
            index = (j-1)*num_replicates_current*num_groups + (k-1)*num_groups + i;
            avg_dist_all(index) = avg_dist;
            time_points_all(index) = j;
            group_all{index} = char(group_names{i});
         end
    end
end

%% Plot Data
% Create a scatter plot of average distance versus time point, colored by treatment group
% Plot Data
% Create a scatter plot of average distance versus time point, colored by treatment group

fig1 = figure;
hold on;

mean_dist = zeros(num_groups, num_time_points);
std_dist = zeros(num_groups, num_time_points);
conversion_factor = 2.33; % 1 pixel = 0.1 micrometers
time_points = unique(time_points_all);

for i = 1:num_groups    % Loop over each treatment group

    % Find the indices of the files corresponding to the current treatment group
    indices = find(strcmp(group_all, group_names{i}));
    if isempty(indices)
        warning('No data points found for treatment group %s', group_names{i});
        continue;
    end

    for j = 1:num_time_points   % Loop over each time point
        % Find the indices of the files corresponding to the current treatment group and time point
        time_indices = find(time_points_all(indices) == j);
        if isempty(time_indices)
            warning('No data points found for treatment group %s at time point %d', group_names{i}, j);
            continue;
        end

        % Extract the average distance for the current treatment group and time point
        dist_time = avg_dist_all(indices(time_indices)) * conversion_factor;

        % Calculate the mean and standard deviation for the current treatment group and time point
        mean_dist(i, j) = mean(dist_time);
        std_dist(i, j) = std(dist_time);
    end

    % Plot the mean and standard deviation for the current treatment group
    errorbar(0:num_time_points-1, mean_dist(i, :), std_dist(i, :), 'o-', 'LineWidth', 1.5);
end

% Add a legend with the name of each treatment group
legend(group_names);

% Add axis labels and a title
xlabel('Time (days)');
xticks(0:num_time_points-1);
ylabel('Average Distance from Center Point (\mum)');
title('Effect of Treatment Groups on Cell Migration');
hold off;

%% Save figure and variables to .mat file
[filepath,name,ext] = fileparts(filename_single);

% Save the variables to a Matlab file
fileName = [name, '_data.mat'];
fileNamePath = fullfile(filepath,fileName);
save(fileNamePath, 'group_all', "group_names","avg_dist_all","time_points_all","index","num_replicates_all","mean_dist","std_dist");

% Add annotations to the saved variables
info.group_all = 'All group names in order';
info.group_names = 'All group names';
info.avg_dist_all = 'All group average distances in order (pixels)';
info.time_points_all = 'All time points in order';
info.index = 'total number of entries';
info.num_replicates_all = 'Number of replicates per time point';
info.mean_dist = 'average distance per treatment group and timepoint';
info.std_dist = 'standard deviation per treatment group and timepoint'

% Save the annotations to the same file
save(fileName, 'info', '-append')

% Export figure 
fileName2 = [name, '_graph.jpg'];
fileNamePath = fullfile(filepath,fileName2);
exportgraphics(fig1,fileNamePath,'BackgroundColor', 'white');
saveas(fig1,fileNamePath,'fig');