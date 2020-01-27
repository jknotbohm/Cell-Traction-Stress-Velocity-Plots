% function plot_cellvel
%PLOT_CELLVEL Plot cell velocities
%
% First run digital image correlation on images of the cells. Since you are
% using image correlation, which correlates a subset, cells must be
% confluent.
%
% This script requires a mat file containing the following:
%   x, y: 2D arrays containing the gridpoints on which the DIC was
%         computed. This can be made using Matlab's meshgrid command.
%         Units: pix.
%   u, v: Incremental displacements of cell image computed by image
%         correlation in horizontal and vertical directions. These are 2D
%         or 3D arrays of size (M, N, P) where M and N are the number of
%         rows and columns, which must match the size of x and y. Variable
%         P corresponds to different time points. If there is only one time
%         point, then the array is 2D (i.e., P=1).
%
% This script also requires a file titled 'ExperimentalSettings.txt'
% containing information described in the readme.
%
% Other required inputs are described in the comments of the section USER
% INPUTS below.
%
% Optional: It can be useful to run batch jobs by running this script as a
% function.
%
% Written by Jacob Notbohm, Univerity of Wisconsin-Madison, 2015-2020


clear;
close all;
clc;


%% --- USER INPUTS ---
% Name of displacment data to load
DICname = 'FIDICc2.mat';
% Name of multipage tif file to plot cells.
%   Set to [] if it matches c2_*.tif
%   Set to 'none' if there is no cell image
cellname = 'c2_island01.tif';
% Name of domain. This is where cells are located. Set to [] if no domain
domainname = 'domain.tif';
% Header of name to save plots. This can contain a directory listing
dirname = 'cell_velocity'; % Name of a folder to put plots in
savenameheader = [dirname,'/t_']; % Header of file name to save
% Set to [] to make figure visible. Set to 1 to make figure invisible.
invisible = [];
% Time between images
time_increment = 10; % min
% Max velocity for color plots
max_vel = 0.25;   % units: um/min
% Typically, the quiver plot, the number of quivers needs to be reduced so
% that they can be seen more clearly. Downsample the number of data points
% for plotting quivers by this factor
qd = 4;


%% --- MAKE PLOTS ---

% Make directory to save plots
if exist(dirname,'dir')==7
    rmdir(dirname,'s');
end
mkdir(dirname);

% Get name of multipagetif file to plot cells
if isempty(cellname)
    cellname = dir('c2_*.tif');
    cellname = cellname(1).name;
end

% Get pixel size from Experimental Settings file
fid = fopen('ExperimentalSettings.txt');
txtcell = cell2mat(textscan(fid,'%f %*[^\n]')); % '%*[^\n]' skips the remainder of each line
pix_size = txtcell(1)*1e6; % Pixel size, um
fclose(fid);
% Load data
load(DICname);
u_cell=u; v_cell=v; x_cell=x; y_cell=y; % Rename variables associated with cell displacements so they aren't overwritten
% Convert from pix to um
x_cell=x_cell*pix_size;     y_cell=y_cell*pix_size;
u_cell=u_cell*pix_size;     v_cell=v_cell*pix_size;
% Convert from displacements to velocities
u_cell = u_cell/time_increment;
v_cell = v_cell/time_increment;

% Number of correlations
K = size(u_cell,3);

% Get center from center of first domain image
if ~isempty(domainname)
    domain1 = imread(domainname,1);
    domain1 = double(domain1); % Convert to double precision
    domain1 = domain1/max(domain1(:)); % Set max value to 1
    % Downsample domain
    % x and y grid points start at w0/2 and end w0/2 before the image ends.
    % First crop off edges so that domain matches start and end points of x
    % and y.
    domain1 = domain1(round(min(y_cell(:))/pix_size):round(max(y_cell(:))/pix_size),round(min(x_cell(:))/pix_size):round(max(x_cell(:))/pix_size));
    domain1 = downsample(domain1,d0); % downsample number of rows
    domain1 = downsample(domain1',d0)'; % downsample number of cols
    % [M, N] = size(x);
    % domain = domain(1:M,1:N); % Correct for slightly larger domain. I should clean this up later.
    % Centroid coordinates
    xc = sum(x(:).*domain1(:)) / sum(domain1(:)); % Units: pix (same units as x and y)
    yc = sum(y(:).*domain1(:)) / sum(domain1(:));
    % Find indices corresponding to nearest x and y coords to xc and yc
    xv = x(1,:);
    yv = y(:,1);
    distx = abs(xc-xv);
    disty = abs(yc-yv);
    [~, xc_idx] = min(distx);
    [~, yc_idx] = min(disty);
    % Convert back to um
    xc = xc_idx*d0*pix_size;
    yc = yc_idx*d0*pix_size;
else
    xc = mean(x(:));
    yc = mean(y(:));
end

cmap = load('map_cold-hot.dat'); % Load hot-cold colormap data

for k=1:K
    hf = make_fig([0.2 0.2 2 1.4]);
    set(hf,'DefaultAxesPosition',[0.03 0.05 .95 .89]);
    if invisible
        set(hf,'visible','off');
    end
    
    % Cell images
    if strcmp(cellname,'none')==0
        % First cell image
        subplot(2,3,1)
        im_k = imread(cellname,k);
        [M, N] = size(im_k);
        imagesc([0 N]*pix_size,[0 M]*pix_size,im_k); colormap(gca,'gray');
        axis xy; axis equal; axis tight; set(gca,'box','off');
        if pix_size == 1
            xlabel('pix'); ylabel('pix');
        else
            xlabel('\mum'); ylabel('\mum');
        end
        title('Reference')
        % Second cell image
        subplot(2,3,4)
        im_k = imread(cellname,k+1);
        [M, N] = size(im_k);
        imagesc([0 N]*pix_size,[0 M]*pix_size,im_k); colormap(gca,'gray');
        axis xy; axis equal; axis tight; set(gca,'box','off');
        xlabel('\mum'); ylabel('\mum');
        title('Current')
    end
    
    % Domain
    if ~isempty(domainname)
        domain = imread(domainname,k);
        domain = double(domain); % Convert to double precision
        domain = domain/max(domain(:)); % Set max value to 1
        domain = logical(domain); % Convert to logical
        % Downsample domain
        % x and y grid points start at w0/2 and end w0/2 before the image ends.
        % First crop off edges so that domain matches start and end points of x
        % and y.
        domain = domain(round(min(y_cell(:))/pix_size):round(max(y_cell(:))/pix_size),round(min(x_cell(:))/pix_size):round(max(x_cell(:))/pix_size));
        domain = downsample(domain,d0); % downsample number of rows
        domain = downsample(domain',d0)'; % downsample number of cols
        % Centroid coordinates
        xc = sum(x(:).*domain(:)) / sum(domain(:)); % Units: pix (same units as x and y)
        yc = sum(y(:).*domain(:)) / sum(domain(:));
        % Find indices corresponding to nearest x and y coords to xc and yc
        xv = x(1,:);
        yv = y(:,1);
        distx = abs(xc-xv);
        disty = abs(yc-yv);
        [~, xc_idx] = min(distx);
        [~, yc_idx] = min(disty);
        % Convert back to um
        xc = xc_idx*d0*pix_size;
        yc = yc_idx*d0*pix_size;
    end
    
    % Cell velocity
    u_cell_k = u_cell(:,:,k);
    v_cell_k = v_cell(:,:,k);
    
    if ~isempty(domainname)
        % Correct for drift by finding mean of velocity outside domain
        SE = strel('disk',5,0);
        domain_dilate = imdilate(domain,SE);
        u_cell_k = u_cell_k - nanmean(u_cell_k(~domain_dilate));
        v_cell_k = v_cell_k - nanmean(v_cell_k(~domain_dilate));
        
        %         % Correct for drift by subtracting off median velocity
        %         u_cell_k = u_cell_k - median(u_cell_k(domain));
        %         v_cell_k = v_cell_k - median(v_cell_k(domain));
        
        % Set values outside domain to zero
        u_cell_k(~domain) = 0;
        v_cell_k(~domain) = 0;
    else
        % Subtract off median of full velocity field. This may or may not
        % work for your images.
        %         u_cell_k = u_cell_k - nanmedian(u_cell_k(:));
        %         v_cell_k = v_cell_k - nanmedian(v_cell_k(:));
    end
    
    % Plot speed
    subplot(2,3,5)
    u_cell_mag = sqrt(u_cell_k.^2 + v_cell_k.^2);
    imagesc([min(x_cell(:)) max(x_cell(:))],[min(y_cell(:)) max(y_cell(:))],u_cell_mag);
    hold on
    x_cell2=downsample(x_cell,qd); x_cell2=downsample(x_cell2',qd)';
    y_cell2=downsample(y_cell,qd); y_cell2=downsample(y_cell2',qd)';
    u_cell_k2=downsample(u_cell_k,qd); u_cell_k2=downsample(u_cell_k2',qd)';
    v_cell_k2=downsample(v_cell_k,qd); v_cell_k2=downsample(v_cell_k2',qd)';
    quiver(x_cell2,y_cell2,u_cell_k2,v_cell_k2, 2.5,... % The scalar number is the relative scaling (length) of the quivers
        'color', [1 1 1]*0.9,'linewidth',1); %You may have to adjust the color of the quivers to show up agaist the colormap
    caxis([0 max_vel]); colormap(gca, 'parula'); colorbar;
    axis xy; axis equal; axis tight; set(gca,'box','off');
    if pix_size == 1
        xlabel('pix'); ylabel('pix');
    else
        xlabel('\mum'); ylabel('\mum');
    end
    title('Cell speed');
    
    % Plot x velocity
    subplot(2,3,2)
    imagesc([min(x_cell(:)) max(x_cell(:))],[min(y_cell(:)) max(y_cell(:))],u_cell_k);
    caxis([-max_vel max_vel]); colormap(gca, cmap); colorbar;
    axis xy; axis equal; axis tight; set(gca,'box','off');
    if pix_size == 1
        xlabel('pix'); ylabel('pix');
    else
        xlabel('\mum'); ylabel('\mum');
    end
    title('X velocity');
    
    % Plot y velocity
    subplot(2,3,3)
    imagesc([min(x_cell(:)) max(x_cell(:))],[min(y_cell(:)) max(y_cell(:))],v_cell_k);
    caxis([-max_vel max_vel]); colormap(gca, cmap); colorbar;
    axis xy; axis equal; axis tight; set(gca,'box','off');
    if pix_size == 1
        xlabel('pix'); ylabel('pix');
    else
        xlabel('\mum'); ylabel('\mum');
    end
    title('Y velocity');
    
    % OPTIONAL: radial and angular compoents of velocity instead of x and y
    %     % Radial velocity
    %     subplot(2,3,2);
    %
    %     % r_grid = sqrt( (x-xc).^2 + (y-yc).^2 ); % Gridpoints
    %     theta = atan2( (y-yc), (x-xc) );
    %
    %     ur = u_cell_k.*cos(theta) + v_cell_k.*sin(theta);
    %
    %     imagesc([min(x_cell(:)) max(x_cell(:))],[min(y_cell(:)) max(y_cell(:))],ur);
    %     caxis([-max_vel max_vel]); colormap(gca, cmap); colorbar;
    %     axis xy; axis equal; axis tight; set(gca,'box','off');
    %     xlabel('\mum','fontsize',11); ylabel('\mum','fontsize',11);
    %     title('Cell radial velocity');
    %
    %     % Angular velocity
    %     subplot(2,3,3)
    %
    %     ut = u_cell_k.*cos(theta+pi/2) + v_cell_k.*sin(theta+pi/2);
    %
    %     imagesc([min(x_cell(:)) max(x_cell(:))],[min(y_cell(:)) max(y_cell(:))],ut);
    %     caxis([-max_vel max_vel]); colormap(gca,cmap); colorbar;
    %     axis xy; axis equal; axis tight; set(gca,'box','off');
    %     xlabel('\mum','fontsize',11); ylabel('\mum','fontsize',11);
    %     title('Cell angular velocity');
    
    % Save
     set(gcf,'PaperPositionMode','auto','InvertHardCopy','off');
    print('-dpng','-r300',[savenameheader,num2str(k,'%0.3d'),'-',num2str(k+1,'%0.3d')]);
    
    close(hf);
    
end