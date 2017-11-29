global gamma beta

%% Files

% For analysing multiple images
% folder= '\\10.0.0.49\pegasus\Pegasusdata\2017_08_30\';
% folder= 'C:\Users\jginer\Dropbox\Postdoc UCLA\Single shot transverse 4D emittance\Emittance Grid Jared\';
% keywords = {'c2_117_58', 'tif'}; % Jared's measurements
% keywords = {'grid400_ S1_1.2_S2_0.345', 'tif'}; %Tokens to look for in the image files.
% keywords = {'grid300_offhotspot_S1_1.18_S2_0.39', 'tif'}; %Tokens to look for in the image files.
% keywords = {'grid300_48mV_10deg_S1-1p21_S2-0p45','tif'};
% keywords = {'grid300_countsyag2_2.3e4_spec1current_1.26_phase_10deg_gl_50mV_ S1_1.25_S2_0.35','tif'};
% keywords = {'grid400_countsyag2_1e4_spec1current_1.26_phase_10deg_gl_50mV_ S1_1.25_S2_0.35','tif'};
% keywords = {'grid400_countsyag2_2.3e4_spec1current_1.26_phase_10deg_gl_50mV_ S1_1.25_S2_0.3375','tif'};
% keywords = {'grid400_countsyag2_3e3_spec1current_1.26_phase_10deg_gl_50mV_ S1_1.25_S2_0.35','tif'};
% keywords = {'grid400_countsyag6_5e4_spec1current_1.26_phase_10deg_gl_50mV_ S1_1.23_S2_0.37','tif'};
% keywords = {'grid2000_countsyag2_2.3e4_spec1current_1.26_phase_10deg_gl_50mV_ S1_1.25_S2_0.3175','tif'};

% file1 = 'grid400_countsyag2_3e3_spec1current_1.26_phase_10deg_gl_50mV_ S1_1.25_S2_0.35 2017 August 30 19_58_51 125.tif'; %Only for a single image
% file1 = 'grid400_countsyag9_5e4_spec1current_1.26_phase_10deg_gl_50mV_ S1_1.21_S2_0.37 2017 August 30 20_10_58 225.tif'; %Only for a single image

%%

folder= '\\10.0.0.49\pegasus\Pegasusdata\2017_11_03\LightField\';
% keywords = {'TEM300_ S1_1.24_S2_0.49_counts_6p7e4','52.tif'}; target = 'TEM';
% keywords = {'TEM300_ S1_1.24_S2_0.49_counts_12e4','tif'}; target = 'TEM';
keywords = {'PepperPotCu300_ S1_1.24_S2_0.49_counts_6p7e4','36.tif'}; target = 'PP';
% keywords = {'PepperPotCu300_ S1_1.24_S2_0.49_counts_12e4','.tif'}; target = 'PP';

% file1 = 'TEM300_ S1_1.24_S2_0.49_counts_12e4 2017 November 03 11_52_29 69.tif'; target = 'TEM'; %Only for a single image

%% Setup parameters (grid+screen)

% pxconv = 12.9e-6; % pixel size [m]  Sept 30th
pxconv = 13.8e-6; % pixel size [m]  Nov 3rd

% gamma = 0.075*50+2.607;
gammabeta = 5.4*(1.71+.05); gamma=sqrt(gammabeta^2+1); % Nov 3rd
beta = sqrt(1-1/gamma^2);

gridSpacing = 85e-6; % [m]
pitch_to_bar_width_ratio = 85/31;
driftLength = 0.85; % [m] drift length from grid to screen


%% Analysis parameters

%maxs = 10;                       % Maximum sigma of the fitted ERFs, in pixels. Larger values will be considered garbage and thrown out.
minpeakdistance = 30;     % In search for peaks, this is the minimum distance allowed between two bars.
maxpeakdistance = 55;
maxtroughheight_x = -0.65;  % This is the minimum height of a peak still considered a 'peak' (TEM grids)
maxtroughheight_y = -0.65;  % This is the minimum height of a peak still considered a 'peak' (TEM grids)
minpeakheight_x = 0.1;  % This is the minimum height of a peak still considered a 'peak' (PepperPot)
minpeakheight_y = 0.1; % This is the minimum height of a peak still considered a 'peak' (PepperPot)
threshint = 0.03; % Threshold of intensity in each beamlet (as a fraction of total intensity)
interval_pc = 60; % interval in [%] of integration along y (x) to fit sigmaxp (sigmayp), centered at midpoint of the bar


%% Crop settings (ellipse)
%     center = [638 505]; % for measurements taken on Aug 30th 2017
%     rx= 290;
%     ry= 295;
    
%     center = [576 534]; % for measurements taken on Aug 18th 2017 and Aug 29th
%     rx= 280;
%     ry= 275;
    
%     center = [535 546]; % for measurements taken on Aug 10th 2017
%     rx= 300;
%     ry= 290;

%     center = [410 580]; % for Jared's measurements
%     rx= 300;
%     ry= 290;

    center = [599 472]; % for measurements taken on Nov 3rd 2017
    rx= 287;
    ry= 289;