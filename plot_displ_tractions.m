%function plot_displ_tractions
%PLOT_DISPL_TRACTIONS Plot displacements and tractions
% 
% First run digital image correlation and compute tractions.
% 
% This script requires a mat file containing the following:
%   x, y: 2D arrays containing the gridpoints on which the DIC was 
%         computed. This can be made using Matlab's meshgrid command. 
%         Units: pix. 
%   u, v: Displacements computed by image correlation in horizontal and 
%         vertical directions. These are 2D or 3D arrays of size (M, N, P) 
%         where M and N are the number of rows and columns, which must 
%         match the size of x and y. Variable P corresponds to different 
%         time points. If there is only one time point, then the array is 
%         2D (i.e., P=1). 
%   tx, ty: Tractions in horizontal and vertical directions. 2D or 3D 
%           arrays of size (M, N, P) that matches the size of u and v.
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
% Set to [] to make figure visible.
invisible = [];
% Name of multipage tif file to plot cells. 
%   Set to [] if 'c2' is in the name
%   Set to 'none' if there is no cell image
cellname = 'c2_island01.tif';
% Name of mat file with displacements and tractions
filename = 'tract_results.mat';
% Name of domain. This is where cells are located. Set to [] if no domain
domainname = [];
% Header of name to save plots.
dirname = 'displ_traction'; % Name of a folder to put plots in
savenameheader = [dirname,'/t_']; % Header of file name to save
% savenameheader = 'displ_traction_'; % Other options are commented
% curdir = pwd;
% [~, folderName, ~] = fileparts(curdir);
% savenameheader = ['../',folderName,'displ_traction_t'];
% Max value of displ and traction (used for color plot limits)
umax = 1;     % um
tmax = 100;      % Pa
% Number of images to plot. Set to empty arry [] to plot all images
num_images = [];

%% --- LOAD DATA ---

% Make directory to plot data
if exist(dirname,'dir')==7
    rmdir(dirname,'s');
end
mkdir(dirname);

% Get name of multipagetif file to plot cells
if isempty(cellname)
    cellname = dir('*c2*.tif');
    cellname = cellname(1).name;
end

% Get pixel size from Experimental Settings file
fid = fopen('ExperimentalSettings.txt');
txtcell = cell2mat(textscan(fid,'%f %*[^\n]')); % '%*[^\n]' skips the remainder of each line
pix_size = txtcell(1)*1e6; % Pixel size, um
fclose(fid);
% Load data
load(filename);
x=x*pix_size; y=y*pix_size;
u=u*pix_size; v=v*pix_size;
% Number of correlations
if isempty(num_images)
    num_images = size(u,3);
end

% Get center from center of first domain image
if ~isempty(domainname)
    domain1 = imread(domainname,1);
    domain1 = double(domain1); % Convert to double precision
    domain1 = domain1/max(domain1(:)); % Set max value to 1
    % Downsample domain
    % x and y grid points start at w0/2 and end w0/2 before the image ends.
    % First crop off edges so that domain matches start and end points of x
    % and y.
    domain1 = domain1(round(min(y(:))/pix_size):round(max(y(:))/pix_size),round(min(x(:))/pix_size):round(max(x(:))/pix_size));
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

% Load hot-cold colormap data
cmap = load('map_cold-hot.dat');

%% --- MAKE PLOTS ---

for k=1:num_images
    hf = make_fig([0.1 0 3 2.4]);
    set(hf,'defaultaxesposition',[0.05 0.08 .92 .86]); % Set axes to nearly fill window
    if invisible
        set(hf,'visible','off');
    end
    
    % Cell image
    if strcmp(cellname, 'none') == 0
        subplot(3,4,1);
        im_k = imread(cellname,k);
        thr = prctile(im_k(:), 99.2); im_k(im_k>thr)=thr;
        thr = prctile(im_k(:), 1); im_k(im_k<thr)=thr;
        [M, N] = size(im_k);
        imagesc([0 N]*pix_size,[0 M]*pix_size,im_k);
        xlabel('\mum'); ylabel('\mum');
        axis xy; axis equal; axis tight; set(gca,'box','off');
        colormap(gca,'gray');
    end
    
    % Displacement Plots
    subplot(3,4,2);
    imagesc([min(x(:)), max(x(:))],[min(y(:)) max(y(:))],u(:,:,k));
    colormap(gca,cmap);
    caxis([-umax umax]); hc = colorbar;
    xlabel('\mum'); ylabel('\mum');
    axis xy; axis equal; axis tight; axis tight; set(gca,'box','off');
    title('U_x')
    
    subplot(3,4,3);
    imagesc([min(x(:)), max(x(:))],[min(y(:)) max(y(:))],v(:,:,k));
    colormap(gca,cmap);
    caxis([-umax umax]); hc = colorbar;
    xlabel('\mum'); ylabel('\mum');
    axis xy; axis equal; axis tight; set(gca,'box','off');
    title('U_y')
    
    subplot(3,4,4)
    imagesc([min(x(:)), max(x(:))],[min(y(:)) max(y(:))],sqrt(u(:,:,k).^2+v(:,:,k).^2));
    colormap(gca,'hot');
    caxis([0 umax]); hc = colorbar;
    xlabel('\mum'); ylabel('\mum');
    axis xy; axis equal; axis tight; set(gca,'box','off');
    title('|U|')
    
    % Traction Plots
    subplot(3,4,6)
    imagesc([min(x(:)), max(x(:))],[min(y(:)) max(y(:))],tx(:,:,k));
    colormap(gca,cmap);
    caxis([-tmax tmax]); hc = colorbar;
    xlabel('\mum'); ylabel('\mum');
    axis xy; axis equal; axis tight; set(gca,'box','off');
    title('t_x')
    
    subplot(3,4,7)
    imagesc([min(x(:)), max(x(:))],[min(y(:)) max(y(:))],ty(:,:,k));
    colormap(gca,cmap);
    caxis([-tmax tmax]); hc = colorbar;
    xlabel('\mum'); ylabel('\mum');
    axis xy; axis equal; axis tight; set(gca,'box','off');
    title('t_y')
    
    subplot(3,4,8)
    imagesc([min(x(:)), max(x(:))],[min(y(:)) max(y(:))],sqrt(tx(:,:,k).^2+ty(:,:,k).^2));
    colormap(gca,'hot');
    caxis([0 tmax]); hc = colorbar;
    xlabel('\mum'); ylabel('\mum');
    axis xy; axis equal; axis tight; set(gca,'box','off');
    title('|t|')
    
%     % OPTIONAL: Quivers to show direction of traction
%     subplot(3,4,5)
%     
%     % plot cell image
%     [M, N] = size(im_k);
%     imagesc([0 N]*pix_size,[0 M]*pix_size,im_k);
%     axis xy;
%     colormap(gca,'gray');
%     hold on;
%     
%     % scale tractions to unit length
%     tx2 = tx(:,:,k)./sqrt(tx(:,:,k).^2 + ty(:,:,k).^2);
%     ty2 = ty(:,:,k)./sqrt(tx(:,:,k).^2 + ty(:,:,k).^2);
%     % downsample data
%     a=6;
%     x2=downsample(x,a); x2=downsample(x2',a)';
%     y2=downsample(y,a); y2=downsample(y2',a)';
%     tx2=downsample(tx2,a); tx2=downsample(tx2',a)';
%     ty2=downsample(ty2,a); ty2=downsample(ty2',a)';
%     % scale to slightly larger than grid spacing
%     scale_factor = 1.1*(x2(1,2)-x2(1,1));
%     % plot quivers
%     quiver(x2,y2,tx2.*scale_factor,ty2.*scale_factor,0,...
%         'showarrowhead','off','linewidth',0.6,'color','g');
%     xlabel('\mum'); ylabel('\mum');
%     axis square; axis tight; set(gca,'box','off');
    
    % OPTIONAL: Convert to radial and hoop components
    theta = atan2(y-yc, x-xc);
    tr = tx(:,:,k).*cos(theta) + ty(:,:,k).*sin(theta);
    ttheta = tx(:,:,k).*cos(theta+pi/2) + ty(:,:,k).*sin(theta+pi/2);
    
    subplot(3,4,10)
    imagesc([min(x(:)), max(x(:))],[min(y(:)) max(y(:))],tr);
    caxis([-tmax tmax]); colormap(gca,cmap); hc = colorbar;
    xlabel('\mum'); ylabel('\mum');
    axis xy; axis equal; axis tight; set(gca,'box','off');
    title('Radial traction, t_r')
    
       
    subplot(3,4,11)
    imagesc([min(x(:)), max(x(:))],[min(y(:)) max(y(:))],ttheta);
    caxis([-tmax tmax]); colormap(gca,cmap); hc = colorbar;
    xlabel('\mum'); ylabel('\mum');
    axis xy; axis equal; axis tight; set(gca,'box','off');
    title('Angular traction, t_\theta')
    
  
    % Save
    set(gcf,'PaperPositionMode','auto');
    savename = [savenameheader,num2str(k,'%3.3d')];
    print(hf,'-dpng','-r200',savename);
    close(hf);
    
end