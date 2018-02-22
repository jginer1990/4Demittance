%%
filename='\\win.desy.de\home\danmarx\My Documents\Project\UCLA\LowEmittance\test_yx_85_31_cor_new.mat';

%% Setup parameters (grid+screen)
pxconv = 15e-6; % pixel size [m]

mask_prop.gridSpacing = 85e-6; % [m]
mask_prop.pitch_to_bar_width_ratio = 85/31;
mask_prop.driftLength = 0.85; % [m] drift length from grid to screen

%% Analysis parameters

%maxs = 10;                       % Maximum sigma of the fitted ERFs, in pixels. Larger values will be considered garbage and thrown out.
analysis.minpeakdistance = 3;     % In search for peaks, this is the minimum distance allowed between two bars.
analysis.maxpeakdistance = 55;
analysis.maxtroughheight_x = -1;  % This is the minimum height of a peak still considered a 'peak' (TEM grids)
analysis.maxtroughheight_y = -1;  % This is the minimum height of a peak still considered a 'peak' (TEM grids)
% analysis.threshint = 0.002; % Threshold of intensity in each beamlet (as a fraction of total intensity)
analysis.threshint = 0.02; % Threshold of intensity in each beamlet (as a fraction of total intensity)
analysis.interval_pc = 60; % interval in [%] of integration along y (x) to fit sigmaxp (sigmayp), centered at midpoint of the bar
analysis.sigma_initguess = 3;