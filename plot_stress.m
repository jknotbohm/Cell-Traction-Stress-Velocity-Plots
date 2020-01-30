% function plot_stress
%PLOT_STRESSE Plot stresses computed by monolayer stress microscopy.
% 
% The input to this script is the mat file output from
% run_stress_calculation.m.
% 
% Additionally, this script requires a mat file containing the following:
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
%   Set to [] if it starts with c2
%   Set to 'none' if there is no cell image
cellname = [];
% Header of name to save plots. This can contain a directory listing
dirname = 'stress_plots';
savenameheader = [dirname,'/t_'];
% curdir = pwd;
% [~, folderName, ~] = fileparts(curdir);
% savenameheader = ['../',folderName,'_stress_t'];
% Name of multipage tif file for domain. Set to [] if monolayer covers
% entire image (ie, domain==1 everywhere).
% domainname = 'domain.tif';
domainname = [];
% Name of mat file with stresses
s_filename = 'stresses.mat';
% Name of mat file with displacements and tractions--needed to get grid
% spacing d0
ut_filename = 'displ_tractions.mat';
% Approx max stress - used for contour limits
smax = 500;  % units: Pa

%% --- LOAD DATA AND PLOT ---

% Make directory to save plots
if exist(dirname,'dir')==7
    rmdir(dirname,'s');
end
mkdir(dirname);

% Get name of multipagetif file to plot cells
if isempty(cellname)
    %     cellname = dir('c2_island*.tif');
    cellname = dir('c2*.tif');
    cellname = cellname(1).name;
end

% Load data
load(s_filename);
load(ut_filename,'d0'); % Get subset spacing

% Get pixel size
fid = fopen('ExperimentalSettings.txt');
txtcell = cell2mat(textscan(fid,'%f %*[^\n]')); % '%*[^\n]' skips the remainder of each line
pix_size = txtcell(1)*1e6; % Pixel size, um
strip_island = txtcell(4); % variable==1 for strip, 2 for island
fclose(fid);

% Convert to um
x=x*pix_size; y=y*pix_size;
% Number of correlations
K = size(Sxx,3);

% Load colormap files
cmap = load('map_cold-hot.dat');
cmap_wrap = load('map_cold-hot-wrap.dat');

for k=1:K
    hf = make_fig([0.1 0.1 3 1.5]);
    if invisible
        set(hf,'visible','off');
    end
    
    % Get k-th domain
    if ~isempty(domainname)
        domain = imread(domainname,k);
        domain = double(domain);
        % Downsample domain
        % x and y grid points start at w0/2 and end w0/2 before the image ends.
        % First crop off edges so that domain matches start and end points of x
        % and y.
        domain = domain(round(min(y(:))/pix_size):round(max(y(:))/pix_size),round(min(x(:))/pix_size):round(max(x(:))/pix_size));
        domain = downsample(domain,d0); % downsample number of rows
        domain = downsample(domain',d0)'; % downsample number of cols
        domain = domain/max(domain(:));
        domain = logical(domain);
    else
        domain = true(size(x));
    end
    
%     % Smooth domain by dilating
%     SE = strel('disk',4,0);
%     domain = imdilate(domain,SE);
    
    % Cell image
    if strcmp(cellname, 'none') == 0
        subplot(2,4,1)
        im_k = imread(cellname,k);
%         % For edge/strip, flip image
%         if strip_island==1
%             domain_col1 = domain(:,1);
%             if mean(domain_col1) < 0.5 % Only flip if domain's 1st col is on right
%                 im_k = fliplr(im_k);
%             end
%         end
        [M, N] = size(im_k);
        thr = prctile(im_k(:),99.9);
        im_k(im_k>thr)=thr;
        thr = prctile(im_k(:),0.1);
        im_k(im_k<thr)=thr;
        imagesc([min(x(:)), max(x(:))],[min(y(:)), max(y(:))],im_k);
        colormap(gca,'gray');
        axis xy; axis equal; axis tight; set(gca,'box','off');
        xlabel('\mum'); ylabel('\mum');
    end
    
    % Stresses
    Sxx_k = Sxx(:,:,k);
    Syy_k = Syy(:,:,k);
    Sxy_k = Sxy(:,:,k);
    
    % Set values outside domain to zero
    Sxx_k(~domain)=0;
    Syy_k(~domain)=0;
    Sxy_k(~domain)=0;
    
    % Principal stresses and directions
    S1 = zeros(size(Sxx_k)); S2=S1; p_angle1=S1; % Preallocate for speed
    [M, N] = size(S1);
    for m=1:M
        for n=1:N
            if ~isnan(Sxx_k(m,n)) && ~isnan(Sxy_k(m,n)) && ~isnan(Syy_k(m,n))
                A = [Sxx_k(m,n) Sxy_k(m,n) ; Sxy_k(m,n) Syy_k(m,n)];
                [v, lambda] = eig(A);
                [S1(m,n), col] = max(diag(lambda));
                S2(m,n) = min(diag(lambda));
                v1 = v(:,col); % Eigenvector corresponding to largest eigenvalue
                %                 p_angle1(m,n) = atan2(v1(2),v1(1));
                p_angle1(m,n) = atan(v1(2)/v1(1));
            else
                S1(m,n)=nan; S2(m,n)=nan; p_angle1(m,n)=nan;
            end
        end
    end
    
    % Mean normal stress
    S_mean = (Sxx_k + Syy_k)/2;
    % Max shear stress
    S_xymax = (S1 - S2)/2;
    
    subplot(2,4,2)
    h1 = pcolor(x,y,S1);
    set(h1,'linestyle','none');
    colormap(gca,cmap);
    caxis([-smax smax]); hc = colorbar;
    axis xy; axis equal; axis tight; set(gca,'box','off');
    xlabel('\mum'); ylabel('\mum');
    title('\sigma_1')
    
    subplot(2,4,3)
    h1 = pcolor(x,y,S2);
    set(h1,'linestyle','none');
    colormap(gca,cmap);
    caxis([-smax smax]); hc = colorbar;
    axis xy; axis equal; axis tight; set(gca,'box','off');
    xlabel('\mum'); ylabel('\mum');
    title('\sigma_2')
    
    subplot(2,4,4)
    h1 = pcolor(x,y,p_angle1*180/pi);
    set(h1,'linestyle','none');
    colormap(gca,cmap_wrap);
    caxis([-90 90]); hc = colorbar;
    axis xy; axis equal; axis tight; set(gca,'box','off');
    xlabel('\mum'); ylabel('\mum');
    title('1st principal angle')
    
    subplot(2,4,6)
    h1 = pcolor(x,y,S_mean);
    set(h1,'linestyle','none');
    caxis([-smax smax]); colormap(gca,cmap); hc = colorbar;
    axis xy; axis equal; axis tight; set(gca,'box','off');
    xlabel('\mum'); ylabel('\mum');
    title('Tension, (\sigma_1+\sigma_2)/2')
    
    subplot(2,4,7)
    h1 = pcolor(x,y,S_xymax);
    set(h1,'linestyle','none');
    colormap(gca,cmap);
    caxis([-smax smax]/2.5); hc = colorbar;
    axis xy; axis equal; axis tight; set(gca,'box','off');
    xlabel('\mum'); ylabel('\mum');
    title('Max shear, (\sigma_1-\sigma_2)/2')
    
%     % OPTIONAL: Ratio of max shear to mean stress
%     subplot(2,4,8)
%     h1 = pcolor(x,y,S_xymax./S_mean);
%     set(h1,'linestyle','none');
%     colormap(gca,'parula');
%     caxis([0 1]); hc = colorbar;
%     axis xy; axis equal; axis tight; set(gca,'box','off');
%     xlabel('\mum'); ylabel('\mum');
%     title('Max shear / Tension')
    
%     % OPTIONAL: Histogram of tension
%     subplot(2,4,5)
%     hist_bins = linspace(0,smax,50);
%     f = hist(S_mean(domain==1),hist_bins);
%     f = f/length(S_mean(domain==1));
%     plot(hist_bins,f*100,'linewidth',2);
%     xlabel('Tension (Pa)')
%     ylabel('%')
%     set(gca,'box','off');
%     xlim([0 smax])
    
    % Save
    set(gcf,'PaperPositionMode','auto','InvertHardCopy','off');
    print('-dpng','-r200',[savenameheader,num2str(k,'%3.3d')]);
    close(hf);
    
end
