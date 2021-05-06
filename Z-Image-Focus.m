%% Image Z Stack Maker
% By R Scales May 2021
%
% This code is not designed as a function. All you mainly need to do is
% click the run button, and change some of the settings in the settings
% section.
%
% If you have any questions, please feel free to contact me on GitHub at
% IonicRob.

clc
fprintf('Started code...\n\n');
clear
close all

% Getting File Locations
[files,path] = uigetfile({'*.png;*.tif'},'MultiSelect','on');
fullfiles = string(fullfile(path,files));
files = string(files);

%% Main Section Part I - Settings
fprintf('Began main section of code...\n\n');
close all
clearvars -except files path fullfiles

% Core Settings
%   Change these to change how it processes the images.
tolerance = 5; % This is the upper limit on how many standard deviations a grid space needs to be to be considered "in focus".
noise_neighbours = 7; % Reco. 7. This is how many nearest neighbours the smoothed image for the grid image selection has to be. If larger the image has larger noise reduction.
col_subdiv = 25; % Reco. 100. This is how many subdivisions are made in the x dimension to form the grids that will be used to build up the final image.
min_subdivision_size = 1; % 1 is the lowest limit, any lower and the pixel grid size would be less than 1 pixel!

% Run Settings
%   Change these to change how it shows how it processes the data.
show_binning = false; % Reco. off. Shows which grids fell under the tolerance set by "tolerance".
show_image_build_up = false; % Reco. off. esp. for large images! Shows the image being built up grid by grid.
assessing_loading_bars_on = false; % Reco. off. Loading bar for it assessing the images. Slows down code.
focusing_loading_bars_on = true; % Reco. off. Loading bar for it processing the final image. Slows down code.
loading_bar_subdivisions = 20; % With the loading bars above, this is the number of divisions you want the loading bars to show. If too large it slows down the code!

%% Main Section Part II - Finding Blurriness

% Useful Variables - Do not alter!
num_of_col_lines = col_subdiv+1; % Number of column lines i.e. the number of x divisons
[aspect_ratio, row_subdiv] = useful_function_1(fullfiles, col_subdiv);
num_of_row_lines = row_subdiv+1; col_subdiv+1; % Number of row lines i.e. the number of y divisons
number_of_regions = col_subdiv*row_subdiv; % This is the total number of grid regions.
num_of_images = length(files);
focus_degree_array = nan(row_subdiv, col_subdiv, num_of_images); % This array is used to store the degree of focus of the different regions in each of the images.

% loading bars progress faster. If too large then speed will be reduced.
check_points = check_points_maker(number_of_regions, loading_bar_subdivisions);

codeStartTime = tic;

for file_num = 1:num_of_images
    fprintf("When file_num = %i, file is %s\n",file_num,files(file_num));
    rgb = imread(fullfiles(file_num));
    try
        grayImage = rgb2gray(rgb);
    catch
        grayImage = im2gray(rgb);
    end
    sdImage = stdfilt(grayImage, true(3));

    binaryImage = sdImage < tolerance; % Whatever works, like 5 or 20 or whatever.

    [rows, columns, ~] = size(sdImage);
    
    [row_positions, row_spacing] = line_point_generator(rows, num_of_row_lines);
    [col_positions, col_spacing] = line_point_generator(columns, num_of_col_lines);
    if min([row_spacing, col_spacing]) <= min_subdivision_size
        [value,index] = min([row_spacing, col_spacing]);
        answer_array = ["Row Spacing", "Column Spacing"];
        message = sprintf('Number of image subdivisions is too much! Must be below %i!\n%s = %0.1f\n%s = %0.1f',min_subdivision_size,answer_array(1),row_spacing,answer_array(2),col_spacing);
        f_error = errordlg(message,'ZStack');
        waitfor(f_error);
        return
    end
    numberOfColorChannels = size(rgb,3);
    fprintf('\tNumber of colour channels in image = %i\n',numberOfColorChannels);
    
    
    if show_binning == true
        imagesc(binaryImage);
        colorbar
        axis image
        hold on;
        for i = 1:length(row_positions)
            line([1, columns], [row_positions(i), row_positions(i)], 'Color', 'r')
        end
        for i = 1:length(col_positions)
            line([col_positions(i), col_positions(i)], [1, rows], 'Color', 'r')
        end
        title(sprintf('Tolerance value =%i',tolerance))
        pause(3);
    end

    if assessing_loading_bars_on == true
        processing_bar = waitbar(0,sprintf('Assessing image num %i... %.3g%%', file_num, 0));
    end
    
    for i = 1:row_subdiv
        for j = 1:col_subdiv
            if assessing_loading_bars_on == true
                tic
            end
            y_range = col_positions(j):col_positions(j+1);
            x_range = row_positions(i):1:row_positions(i+1);
            blurAreaFraction = nnz(binaryImage(x_range,y_range)) / numel(binaryImage(x_range,y_range));
            focus_degree_array(i,j,file_num) = blurAreaFraction;
            current_point = ((i-1)*col_subdiv)+j;
            YN = plot_current_progress_yn(current_point, check_points);
            if assessing_loading_bars_on == true && YN == true
                percentage = current_point/(number_of_regions);
                waitbar(percentage,processing_bar,sprintf('Assessing image num %i...\n%.3g%%', file_num, percentage*100));
            end
        end
    end
    if assessing_loading_bars_on == true
        close(processing_bar)
    end
    clear rgb grayImage current_point YN sdImage binaryImage
end

fprintf('Finished evaluating images...\n\tTook %.3g seconds to complete!\n\n',toc(codeStartTime));
clear blurAreaFraction
%% %% Main Section Part III - Filtering
fprintf('Beginning image selection and noise reduction section...\n');
codeStartTime = tic;

[M,I] = min(focus_degree_array,[],3);
I_old = I;
I = medfilt2(I,[noise_neighbours,noise_neighbours],'symmetric');
clear focus_degree_array
%% %% Main Section Part IV - Focusing

fprintf('Beginning image focusing section...\n');
close all

final_image = zeros(rows,columns,numberOfColorChannels);
fig_final_image = figure('Name','Final Image');

if focusing_loading_bars_on == true
    loading_bar = waitbar(0,sprintf('Focusing image... %.3g%%',0));
end

images = cell([length(fullfiles),1]);
for u = 1:length(fullfiles)
    images{u} = imread(fullfiles(u));
end


for i = 1:row_subdiv
    for j = 1:col_subdiv
        tic
        y_range = col_positions(j):col_positions(j+1);
        x_range = row_positions(i):1:row_positions(i+1);
        best_image_index = I(i,j);
        best_image = images{best_image_index}; % imread(fullfiles(best_image_index));
        final_image(x_range,y_range,:) = best_image(x_range,y_range,:);
        if show_image_build_up == true
            figure(fig_final_image);
            imshow(uint8(final_image));
            axis image
        end
        current_point = ((i-1)*col_subdiv)+j;
        YN = plot_current_progress_yn(current_point, check_points);
        if focusing_loading_bars_on == true && YN == true
            percentage = current_point/(number_of_regions);
            waitbar(percentage,loading_bar,sprintf('Focusing image...\n%.3g%%', percentage*100));
        end

        clear best_image
    end
end

%% %% Main Section Part V - Final Section

plot_mapping_map = true;
if plot_mapping_map == true
    [X,Y] = meshgrid(col_positions(1:length(col_positions)-1),row_positions:length(row_positions)-1);
    fig_mapping_map = figure('Name','Mapping Map','WindowState','Maximized');
    subplot(1,2,2);
    h1 = pcolor(X,Y,flipud(I));
    set(h1, 'EdgeColor', 'none');
    ax1 = gca;
    ax1.PlotBoxAspectRatio = [1,aspect_ratio,1];
    title(sprintf('Noise Reduction (NN = %i)', noise_neighbours));
    subplot(1,2,1);
    h2 = pcolor(X,Y,flipud(I_old));
    set(h2, 'EdgeColor', 'none');
    ax2 = gca;
    ax2.PlotBoxAspectRatio = [1,aspect_ratio,1];
    title(sprintf('Original Image Selection Map'));
end

if focusing_loading_bars_on == true
    close(loading_bar);
end

if show_image_build_up == false
    figure(fig_final_image);
    imshow(uint8(final_image));
    axis image
end

save_name = sprintf('focused_image tol_%i col_%i row_%i noise_%i.png', tolerance, col_subdiv, row_subdiv, noise_neighbours);
imwrite(uint8(final_image),fullfile(path,save_name),'png');
save_name = sprintf('image_map tol_%i col_%i row_%i noise_%i.png', tolerance, col_subdiv, row_subdiv, noise_neighbours);
saveas(fig_mapping_map,fullfile(path,save_name));

codeEndTime = toc(codeStartTime);
fprintf('Finished focusing images...\n\tTook %.3g seconds to complete!\n\n',toc(codeStartTime));
%% Internal Functions

function [aspect_ratio, row_subdiv] = useful_function_1(fullfiles, col_subdiv)
    rgb_test = imread(fullfiles(1));
    aspect_ratio = size(rgb_test,1)/size(rgb_test,2);
    row_subdiv = round(aspect_ratio*col_subdiv,0); 
end

function check_points = check_points_maker(total_number_of_points, number_of_divisions)
    exact_check_points = linspace(0, total_number_of_points, number_of_divisions);
    check_points =  floor(exact_check_points);
end

function [line_positions, line_spacing] = line_point_generator(size, num_of_lines)
    raw_positions = linspace(1,size,num_of_lines);
    line_positions = floor(linspace(1,size,num_of_lines));
    line_spacing = raw_positions(2)-raw_positions(1);
end

function YN = plot_current_progress_yn(current_point, check_points)
    YN = any(check_points(:) == current_point); % (current_point == check_points);
end
