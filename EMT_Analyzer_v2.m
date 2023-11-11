%% Preamble
clearvars;
close all;

% Get the number of time points from the user
prompt = {'Enter the number of time points:'};
dlgtitle = 'Time points';
dims = [1 35];
definput = {'1'};
answerTP = inputdlg(prompt,dlgtitle,dims,definput);
num_time_points = str2double(answerTP{1});

% Loop over each time point and ask the user to input the number of replicates
for j = 1:num_time_points
    num_replicates = inputdlg(sprintf('Enter the number of replicates for time point %d', j), sprintf('Replicates for Time Point %d', j), 1);
    num_replicates_all(j) = str2double(num_replicates{1});
end

% Prompt the user to input the number and names of treatment groups
prompt = {'Enter the number of treatment groups:'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'1'};
answerTG = inputdlg(prompt, dlgtitle, dims, definput);
num_groups = str2double(answerTG{1});
group_names = cell(num_groups,1);


numberFiles = (num_time_points+length(num_replicates_all))*num_groups;

% Initialize vectors to store the data for each group
avg_dist_all = zeros(numberFiles, 1);
time_points_all = zeros(numberFiles, 1);
group_all = cell(numberFiles, 1);

%% Main Script


for i = 1:num_groups    % Loop over each treatment group

    prompt = sprintf('Enter the name of treatment group %d:', i);
    group_names{i} = char(inputdlg(prompt, dlgtitle, dims));
     
    for j = 1:num_time_points   % Loop over each time point and replicate
        
%             % Prompt the user to input the number of replicates for the current time point
%             num_replicates_current = inputdlg(sprintf('Enter the number of replicates for time point %d:', j));
%             % Convert the input to a numeric value
%              num_replicates_current = str2double(num_replicates_current{1});
        for k = 1:num_replicates_all(j)
           
        
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
            index = (j-1)*num_replicates_all(j)*num_groups + (k-1)*num_groups + i;
            avg_dist_all(index) = avg_dist;
            time_points_all(index) = j;
            group_all{index} = char(group_names{i});
         end
    end
end

%% Plot Data
% Create a scatter plot of average distance versus time point, colored by treatment group

fig1 = figure;
hold on;

mean_dist = zeros(1, length(num_replicates_all));
std_dist = zeros(1, length(num_replicates_all));
conversion_factor = 2.33; % 1 pixel = 0.1 micrometers
time_points = unique(time_points_group);

for i = 1:length(group_names)
    % Find the indices of the files corresponding to the current treatment group
    indices = zeros(index,1);
    for j = 1:index
        indices(j,1) = strcmp(group_all{j,1}, group_names{i});
        % Check if any data points were found for the current treatment group
    end
    indices = find(indices);
    if ~any(indices)
            warning('No data points found for treatment group %s', group_names{i});
            continue;
    end
    % Extract the average distance and time point for the current treatment group
    dist_group = avg_dist_all(indices)* conversion_factor;
    time_points_group = time_points_all(indices);

    % calculate the mean and standard deviation for each time point

    start_idx = 1;
    for s = 1:length(num_replicates_all)
        end_idx = start_idx + num_replicates_all(s) - 1;
        mean_dist(i,s) = mean(dist_group(start_idx:end_idx));
        std_dist(i,s) = std(dist_group(start_idx:end_idx));
        start_idx = end_idx + 1;
    end

    % plot the mean and standard deviation for each time point
    errorbar(0:length(time_points)-1, mean_dist(i,:), std_dist(i,:),'o-', 'LineWidth', 1.5);

    % Plot the data for the current treatment group
%     scatter(time_points_group, avg_dist_group, 'filled');
    % Add a legend with the name of each treatment group
    legendname{i} = char(group_names{i,1});
    
end
legend(legendname);

% Add axis labels and a title
xlabel('Time (days)');
xticks(0:1:length(time_points)-1);
ylabel('Average Distance from Center Point (\mum)');
title('Effect of Treatment Groups on Cell Migration');
hold off;

%% Save figure and variables to .mat file
[filepath,name,ext] = fileparts(filename_single);

% Save the variables to a Matlab file
fileName = [name, '_data.mat'];
fileNamePath = fullfill(filepath,fileName);
save(fileNamePath, 'group_all', "group_names","avg_dist_all","time_points_all","index","num_replicates_all");

% Add annotations to the saved variables
info.group_all = 'All group names in order';
info.group_names = 'All group names';
info.avg_dist_all = 'All group average distances in order (pixels)';
info.time_points_all = 'All time points in order';
info.index = 'total number of entries';
info.num_replicates_all = 'Number of replicates per time point';

% Save the annotations to the same file
save(fileName, 'info', '-append')

% Export figure 
fileName2 = [name, '_graph.jpg'];
fileNamePath = fullfill(filepath,fileName2);
exportgraphics(fig1,fileNamePath,'BackgroundColor', 'white');