% function compute_cell_trajectories
% compute_cell_trajectories
%
% Compute pathlines for cell trajectories based on cell displacements
% calculated from image correlation.
% 
% Required inpus are described in the next section USER INPUTS.
% Additionaly, a file called ExperimentalSettings.txt is required, which is
% described in the repository used for computing cell tractions and
% stresses.
% 
% Systematic drift will create errors in the computed trajectories. See the
% lines in the middle of the for loop for options to correct for drift. Use
% whichever is best for your data or write your own method of drift
% correction.
% 
% This script requires inpaint_nans.m, which is available from the Matlab
% file repository.
% 
% Optional: It can be useful to run batch jobs by running this script as a
% function.
%
% This script adapted from code written by Chan Young Park, Harvard School
% of Public Health, 2012
% 
% Written by Jacob Notbohm, University of Wisconsin-Madison, 2015-2020

clear;
close all;
clc;


%% --- USER INPUTS ---

% Name of displacment data to load. Required variables in the mat file:
%   x, y: 2D arrays of grid points corresponding to the center of each
%         subset
%   u, v: 3D arrays of cell displacement at each grid point over time.
DICname = 'FIDICc2.mat';
% Name of file containing phase contrast images used for correlation. Use a
% multipage tif format to save multiple images over time as one file. 
% Use * for wild, and the code will open the first image matching the given 
% format
cellname = 'c1*.tif';
% It is common for the resolution of the image correlation to be finer than
% the cell size. The trajectories should correlate to a cell, so downsample
% the data to match the number of grid points to the number of cells. Use a
% value of 1 to avoid downsampling. Must be positive integer.
fd = 2;
% Name of multipage tif file with domains. Use [] if domain is entire image
domainname = 'domain.tif';
% State whether geometry is an island. Set this to 1 if yes. For island
% geometry, code will set origin to be the center of the island
isisland = 1;
% Optional: choose starting and ending time points for analysis. Enter
% empty array [] to begin at first time point and/or end at last time point
tstart = 11;
tend = [];
% Name to save mat file of data
savename = 'cell_trajectories.mat';
% Name to save plot of trajectories
plotname = 'Trajectories';


%% --- GET TRAJECTORIES ---

% Get pixel size
fid = fopen('ExperimentalSettings.txt');
txtcell = cell2mat(textscan(fid,'%f %*[^\n]')); % '%*[^\n]' skips the remainder of each line
pix_size = txtcell(1)*1e6; % Pixel size, um

% Get cell name. This is needed if a * is used in cell name
F = dir(cellname);
cellname = F(1).name;

% Load data
load(DICname);
% Rename variables associated with cell displacements so they aren't
% confused with substrate displacements
u_cell=u; v_cell=v; x_cell=x; y_cell=y;
% Convert from pix to um
x_cell=x_cell*pix_size;     y_cell=y_cell*pix_size;
u_cell=u_cell*pix_size;     v_cell=v_cell*pix_size;

% Check if domain is entire image
if isempty(domainname)
%     cellname = dir('c2_*.tif');
%     cellname = dir('pos*.tif');
%     cellname = cellname(1).name;
    info = imfinfo(cellname);
    domain1 = ones(info(1).Height,info(1).Width);
else % Otherwise load domain
    domain1 = imread(domainname,1);
    domain1 = double(domain1); % Convert to double precision
    domain1 = domain1/max(domain1(:)); % Set max value to 1
end

% Downsample domain
% x and y grid points don't necessarily start at 1.
% First crop off edges so that domain matches start and end points of x
% and y.
domain1 = domain1(min(y(:)):max(y(:)), min(x(:)):max(x(:))); % Match starting row and col to y and x
domain1 = downsample(domain1,d0); % downsample number of rows
domain1 = downsample(domain1',d0)'; % downsample number of cols

% If geometry is island, shift origin so it's at center of domain
if isisland == 1
    % Get center from center of first domain image
    % Centroid coordinates
    xc = sum(x(:).*domain1(:)) / sum(domain1(:)); % Units: pix (same units as x and y)
    yc = sum(y(:).*domain1(:)) / sum(domain1(:));
    % Convert to um
    xc = xc*pix_size;
    yc = yc*pix_size;
    
    % Shift data so origin is at center of domain
    x_cell = x_cell - xc; % Units: um
    y_cell = y_cell - yc;
end

% Get initial grid points for trajectory computation
x_cell2=downsample(x_cell,fd); x_cell2=downsample(x_cell2',fd)';
y_cell2=downsample(y_cell,fd); y_cell2=downsample(y_cell2',fd)';


% Get domain for first time point
% For an island, find the largest connected region
if isisland == 1
    BNDRY = bwboundaries(domain1); % This should be 1 cell with boundary coordinates
    LB = cellfun(@length,BNDRY);
    [~, idx] = max(LB);
    BNDRY = BNDRY{idx};
    x_bndry = BNDRY(:,2); % x coordinates are the columns
    y_bndry = BNDRY(:,1); % y coordinates are the rows
    % Convert from d0-spaced pix to um
    x_bndry=x_bndry*pix_size*d0; y_bndry=y_bndry*pix_size*d0;
    % Shift so origin is at zero
    x_bndry = x_bndry-xc; % Units: um
    y_bndry = y_bndry-yc;
    
    % Get indices of points inside region given by (x_bndry,y_bndry)
    IDX = inpolygon(x_cell2,y_cell2,x_bndry,y_bndry);
    
% For edge geometry, just use the given domain    
else
    IDX = logical(domain1d);
end


% Determine starting time point for analysis. If first images are all
% zeros, then starting timepoint should be 2.
if isempty(tstart)
    u1 = u_cell(:,:,1);
    v1 = v_cell(:,:,1);
    if sum(u1(:))==0 && sum(v1(:))==0
        tstart=2;
    else
        tstart=1;
    end
end

% Determine final time point for analysis
if isempty(tend)
    tend = size(u_cell,3)-1;
end

% Number of trajectories is given by number of nonzero elements in IDX
num_traj = sum(double(IDX(:)));

% Preallocate arrays for trajectories. These arrays start at each 
% downsampled  grid point within the domain. Rows are different traj and 
% columns are different times.
traj_x = nan*zeros(num_traj,tstart-tend+1);
traj_y = nan*zeros(num_traj,tstart-tend+1);

% First set of trajectories are given by points inside domain (ie, where IDX==1)
traj_x(:,1) = x_cell2(IDX);
traj_y(:,1) = y_cell2(IDX);

% Scan through remaining time points adding to the path coordinates
for k= tstart:tend
    if ~isempty(domainname)
        % Get domain. There's an ambiguity here, because k is correlation
        % number, and the correlation is between the k-th and (k+1)-th
        % images if tstart==1 or the (k-1)-th and k-th images if tstart==2.
        % Here, I'm choosing the starting image as k-tstart+1
        domain = imread(domainname, k-tstart+1);
        
%         % Uncomment this and comment above if domain is only for first time point
%         domain = imread(domainname);
        
        domain = double(domain);
        domain = domain/max(domain(:));
        domain = logical(domain);
        % Downsample domain
        % x and y grid points start at w0/2 and end w0/2 before the image ends.
        % First crop off edges so that domain matches start and end points of x
        % and y.
        domain = domain(min(y(:)):max(y(:)), min(x(:)):max(x(:)));
        domain = downsample(domain,d0); % downsample number of rows
        domain = downsample(domain',d0)'; % downsample number of cols
    else
        domain = true(size(x));
    end
    
    % --- k-th cell displacements ---
    u_cell_k = u_cell(:,:,k); % units: um
    v_cell_k = v_cell(:,:,k);
    
    % --- [OPTIONAL] Remove displacements that are too large ---
    % Set a threshold based on size of image
    xymax = sqrt(max(x_cell(:))*max(y_cell(:))); % Typical image size. Units: um
    thr = 0.2*xymax; % Choose threshold to be some fraction of image size
    idx = abs(u_cell_k)>thr | abs(v_cell_k)>thr;
    u_cell_k(idx) = nan;
    v_cell_k(idx) = nan;
    u_cell_k = inpaint_nans(u_cell_k);
    v_cell_k = inpaint_nans(v_cell_k);
    
    % --- Correct for drift. Choose one of various options ---
    
%     % Correct for drift by subtracting off mean displacement of the slowest
%     % third of the cells
%     idx = ~isnan(u_cell_k) & ~isnan(v_cell_k);
%     u_tmp = u_cell_k(idx);
%     v_tmp = v_cell_k(idx);
%     umag_tmp = sqrt(u_tmp.^2+v_tmp.^2);
%     umag30 = prctile(umag_tmp,30);
%     idx = umag_tmp < umag30;
%     uave = mean( u_tmp(idx) );
%     vave = mean( v_tmp(idx) );
    
%     % Correct for drift by subtracting off medians
%     uave = nanmedian(u_cell_k(:));
%     vave = nanmedian(v_cell_k(:));
    
    % Correct for drift using domain
    SE = strel('disk',5,0);
    domain_dilate = imdilate(domain,SE);
    uave = nanmean(u_cell_k(~domain_dilate));
    vave = nanmean(v_cell_k(~domain_dilate));
    
    u_cell_k = u_cell_k-uave;
    v_cell_k = v_cell_k-vave;
    
    % Set data outside domain to nan
    u_cell_k(~domain)=nan;
    v_cell_k(~domain)=nan;
    
    % Interpolate k-th displacements to gridpoints of previous timepoint
    displ_x = griddata(x_cell,y_cell,u_cell_k,traj_x(:,k-tstart+1),traj_y(:,k-tstart+1));
    displ_y = griddata(x_cell,y_cell,v_cell_k,traj_x(:,k-tstart+1),traj_y(:,k-tstart+1));
        
    % Add to trajectory arrays
    traj_x(:,k-tstart+2) = traj_x(:,k-tstart+1) + displ_x; % Units: um
    traj_y(:,k-tstart+2) = traj_y(:,k-tstart+1) + displ_y;
end

% Save data
save(savename,'traj_x','traj_y'); % Units: um

%% --- PLOT TRAJECTORIES ---

load(savename);

hf = make_fig([0.5 0.5 1 1]);
% Set axes to nearly fill window
set(hf,'DefaultAxesPosition',[0.1 0.1 .86 .86]);
hold on

% % Downsample for plotting
% traj_x = downsample(traj_x,10);
% traj_y = downsample(traj_y,10);

% % Remove some time points
% traj_x = traj_x(:,1:60);
% traj_y = traj_y(:,1:60);

% Choose color by plotting in a loop
% for k=1:num_traj
%     plot(traj_x(k,:),traj_y(k,:),'b')
% end

% Let Matlab automatically choose color
plot(traj_x',traj_y')
xlabel('\mum','fontsize',11);
ylabel('\mum','fontsize',11);
axis equal;
set(gca,'box','on','fontsize',11);
% Manually set axes
% axis([-800 800 -800 800]);
% axis([-500 500 -500 500]);
% set(gca,'xtick',[-500:250:500],'ytick',[-500:250:500]);

set(hf,'Paperpositionmode','auto');
print('-dpng','-r300',plotname);

% Save in different directory
% curdir = pwd;
% print('-dpng','-r300',['../Trajectories/',curdir(end-4:end)]);





