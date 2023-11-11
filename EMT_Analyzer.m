% % Prompt the user to select the input images
% [filename, pathname] = uigetfile({'*.jpg;*.png;*.tif', 'Image files (*.jpg, *.png, *.tif)'}, ...
%     'Select input images', 'MultiSelect', 'on');
% 
% % Check if the user canceled the file selection dialog
% if isequal(filename,0) || isequal(pathname,0)
%     disp('User canceled the file selection dialog');
%     return;
% end
% 
% % Make sure filename is a cell array
% if ~iscell(filename)
%     filename = {filename};
% end

% Get the number of time points from the user
prompt = {'Enter the number of time points:'};
dlgtitle = 'Time points';
dims = [1 35];
definput = {'1'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
num_time_points = str2double(answer{1});

% Get the number of replicates for each time point from the user
prompt = {'Enter the number of replicates for each time point:'};
dlgtitle = 'Replicates';
dims = [1 45];
definput = repmat({'1'}, 1, num_time_points);
answer = inputdlg(prompt,dlgtitle,dims,definput);
num_replicates = str2double(answer);

% Initialize a vector to store the average distance for each file
avg_dist_all = zeros(length(filename), 1);
% Initialize a vector to store the time points for each file
time_points_all = zeros(length(filename), 1);

% Loop over each time point and replicate
for j = 1:num_time_points
    for k = 1:num_replicates
        % Prompt the user to select the input image
        [filename_single, pathname_single] = uigetfile({'*.jpg;*.png;*.tif', 'Image files (*.jpg, *.png, *.tif)'}, ...
            sprintf('Select input image for time point %d, replicate %d', j, k), 'MultiSelect', 'off');
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

        % Compute the distance from the center point to every point on the outline
        distances = sqrt((xy(:,1)-x).^2 + (xy(:,2)-y).^2);

        % Compute the average distance
        avg_dist = mean(distances);

        % Store the average distance and time point for the current file
        index = (j-1)*num_replicates + k;
        avg_dist_all(index) = avg_dist;
        time_points_all(index) = j;
    end
end

% Loop over each time point
avg_dist_by_time_point = zeros(num_time_points, 1);
std_by_time_point = zeros(num_time_points, 1);
file_counter = 1;
for j = 1:num_time_points
    % Compute the average distance and standard deviation for the current time point
    avg_dist_replicate = zeros(num_replicates(j), 1);
    for k = 1:num_replicates(j)
        avg_dist_replicate(k) = calculate_avg_dist(avg_dist_all{file_counter}, center_point);
        file_counter = file_counter + 1;
    end
    avg_dist_by_time_point(j) = mean(avg_dist_replicate);
    std_by_time_point(j) = std(avg_dist_replicate);
end

% Plot the results using errorbar
errorbar(1:num_time_points, avg_dist_by_time_point, std_by_time_point, '-o');

% Label the x-axis with the time points
xticklabels(arrayfun(@(x) sprintf('Time point %d', x), 1:num_time_points, 'UniformOutput', false));

% Set the title and axis labels
title('Average distance from center to outline by time point');
xlabel('Time point');
ylabel('Average distance');
