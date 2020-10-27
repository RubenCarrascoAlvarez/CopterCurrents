function STCFIT_multi_K = generate_STCFIT_from_IMG_SEQ_multi_K(IMG_SEQ, sq_size_m, sq_dist_m,...
         mask_2D,nan_percentage_thr,water_depth_mask_2D,Ux_limits_FG, Uy_limits_FG,U_FG_res,...
         w_width_FG,U_SG_res,w_width_SG, waveLength_limits_m, wavePeriod_limits_sec,...
         K_steps,max_current_shear_expected)     
% Generate STCFIT, the structure which define the windows to run the
% current fit for multiK fit.
%
% fit, and the paremeters.
%   Input:
%       IMG_SEQ: Georferenceced Image sequence structure
%       sq_size_m: square size in meter for the current fit windows.
%       sq_dist_m: square distance in meter between current fit windows.
%       mask_2D:  2d logic mask defining the area useable for the fit [usable == 1]
%       nan_percentage_thr: maximum percentage of NaN (or not usable data) data to enable fit window
%       Ux_limits_FG: [min max] Ux current interval for the first guess [m/s]
%       Uy_limits_FG: [min max] Uy current interval for the first guess [m/s]
%       U_FG_res: resolution in [m/s] for the first guess vector
%       w_width_FG: filter with in rad/s for first guess
%       U_SG_res: resolution in [m/s] for the second guess vector
%       w_width_SG: filter with in rad/s for second guess
%       waveLength_limits_m: wave-lenght limits [min max] in meters, for the dispersion relation fit.
%       wavePeriod_limits_sec: wave period limits [min max] in seconds, for the dispersion relation fit.
%       K_steps: 2xN K steps for depth fit. 
%       max_current_shear_expected: max current shear expected in m/s
%                                   used for the second guess
%
%   Output:
%       STCFIT: Structure defining the current fit input parameters
%
%         STCFIT.Generic.gridX: 2D X grid (Horizontal) in meters
%         STCFIT.Generic.gridY: 2D Y grid (Vertical) in meters 
%         STCFIT.Generic.image: image example
%         STCFIT.Generic.Longitude: Longitude [deg]
%         STCFIT.Generic.Latitude:  Latitude [deg]
%         STCFIT.Generic.heading:   Angle to North [deg]
%         STCFIT.Generic.altitude:  altitude in meters
%         STCFIT.Generic.time_stamp: time stamps in datenum format
%         STCFIT.Generic.usable_data_mask_2D: usable points for fit.
%         STCFIT.Generic.water_depth_mask_2D: water depth.
%
%         STCFIT.Windows.N_fit_windows: number of windows used for the current fit.
%         STCFIT.Windows.w_corners_dim1: 3xN_fit_windows
%         STCFIT.Windows.w_corners_dim2: 3xN_fit_windows
%         STCFIT.Windows.sq_size_m: square size in meter for the current fit windows.
%         STCFIT.Windows.sq_dist_m: square distance in meter between current fit windows.
%         STCFIT.Windows.average_depth_win: average depth in meters
% 
%         STCFIT.fit_param.Ux_FG_2D: MxN Ux first guess matrix
%         STCFIT.fit_param.Uy_FG_2D: MxN Uy first guess matrix
%         STCFIT.fit_param.Ux_SG_2D: MxN Ux second guess matrix (offset matrix)
%         STCFIT.fit_param.Uy_SG_2D: MxN Uy second guess matrix (offset matrix) 
%         STCFIT.fit_param.w_width_FG: first guess filter width in w [rad/s] 
%         STCFIT.fit_param.w_width_SG: second guess filter width in w [rad/s]
%         STCFIT.fit_param.waveLength_limits_m: [min max] wavelength to use [meters]
%         STCFIT.fit_param.wavePeriod_limits_sec: [min max] wave Period to use [seconds]
%         STCFIT.fit_param.K_limits: [min max] wave number to use [rad/m]
%         STCFIT.fit_param.W_limits: [min max] wave frequency to use [rad/sec]
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% (C) 2019, Ruben Carrasco <ruben.carrasco@hzg.de>
%
% This file is part of CopterCurrents
%
% CopterCurrents has been developed by Department of Radar Hydrography, at
% Helmholtz-Zentrum Geesthacht Centre for Materials and Coastal Research 
% (Germany), based in the work exposed in: 
% 
% M. Streßer, R. Carrasco and J. Horstmann, "Video-Based Estimation 
% of Surface Currents Using a Low-Cost Quadcopter," in IEEE Geoscience 
% and Remote Sensing Letters, vol. 14, no. 11, pp. 2027-2031, Nov. 2017.
% doi: 10.1109/LGRS.2017.2749120
%
%
% CopterCurrents is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% CopterCurrents is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with CopterCurrents.  If not, see <http://www.gnu.org/licenses/>.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


size_IMG_SEQ = size(IMG_SEQ.IMG);

if exist('sq_size_m','var') ==0 || isempty(sq_size_m)
    warning('generate_STCFIT_from_IMG_SEQ:sq_size_m not defined');
    % set sq_size_m to fit 10 squares in dim1 (x)
    sq_size_m = (max(IMG_SEQ.gridX(:)) - min(IMG_SEQ.gridX(:))) / (10.2);
end

if exist('sq_dist_m','var') ==0 || isempty(sq_dist_m)
    % set square distance = to square size
    sq_dist_m = sq_size_m;
end

if exist('mask_2D','var') ==0 || isempty(mask_2D)
    % set mask_2D to 1
    mask_2D = true(size_IMG_SEQ(1),size_IMG_SEQ(2));
end

if exist('nan_percentage_thr','var') ==0 || isempty(nan_percentage_thr)
    % set nan_percentage_thr to 5%
    nan_percentage_thr = 5;
end

if exist('water_depth_mask_2D','var') ==0 || isempty(water_depth_mask_2D)
    water_depth_mask_2D = 10;
end

if isscalar(water_depth_mask_2D) % if the mask a unique 1 value, use this value for all the 2D space
   water_depth_mask_2D = repmat(water_depth_mask_2D,[size(IMG_SEQ.IMG,1) size(IMG_SEQ.IMG,2)]);
end

if exist('Ux_limits_FG','var') ==0 || isempty(Ux_limits_FG)
    Ux_limits_FG = [-1.5 1.5];
end

if exist('Uy_limits_FG','var') ==0 || isempty(Uy_limits_FG)
    Uy_limits_FG = [-1.5 1.5];
end

if exist('U_FG_res','var') ==0 || isempty(U_FG_res)
    U_FG_res = 0.1;
end

if exist('w_width_FG','var') ==0 || isempty(w_width_FG)
    w_width_FG = 1;
end

if exist('U_SG_res','var') ==0 || isempty(U_SG_res)
    U_SG_res = 0.025;
end

if exist('w_width_SG','var') ==0 || isempty(w_width_SG)
    w_width_SG = w_width_FG/2;
end

if exist('waveLength_limits_m','var') ==0 || isempty(waveLength_limits_m)
    % set waveLength_limits_m acording to sq_size_m
    pixel_size = sqrt(IMG_SEQ.dx.^2 + IMG_SEQ.dy.^2);
    waveLength_limits_m = [ 4*pixel_size sq_size_m/4];
end

if exist('wavePeriod_limits_sec','var') ==0 || isempty(wavePeriod_limits_sec)
    % set wavePeriod_limits_sec
    % wavePeriod_limits_sec = [ 4*IMG_SEQ.dt  IMG_SEQ.dt*size_IMG_SEQ(3)/4];
    wavePeriod_limits_sec = [ IMG_SEQ.dt  IMG_SEQ.dt*size_IMG_SEQ(3)];
end

if exist('max_current_shear_expected','var') ==0 || isempty(max_current_shear_expected)
    max_current_shear_expected = 0.25;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% obtain windows fit corners according to sq_size_m & sq_dist_m   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sq_size_pix = sq_size_m/IMG_SEQ.dx;
sq_dist_pix = sq_dist_m/IMG_SEQ.dx;

% creating starting and end point windows in dim1
% w_corners_dim1 = center_dim1 ini_dim1 end_dim1

% offset_center_dim1 = size_IMG_SEQ(1) - ((1+floor(size_IMG_SEQ(1)/sq_size_pix))*sq_size_pix);

aux = 1:sq_size_pix:size_IMG_SEQ(1);
offset_center_dim1 = (size_IMG_SEQ(1)-aux(end))/2;
center_dim1 = offset_center_dim1 -(sq_dist_pix/2):sq_dist_pix:size_IMG_SEQ(1); 
ini_dim1 = round(center_dim1 - (sq_size_pix/2));
end_dim1 = round(center_dim1 + (sq_size_pix/2)-1);
center_dim1 = round(center_dim1);

% delete values out of the IMG_SEQ bounds (dim1)
ind2del = ini_dim1<1 | end_dim1>size_IMG_SEQ(1);
center_dim1(ind2del) = [];
ini_dim1(ind2del) = [];
end_dim1(ind2del) = [];


% creating starting and end point windows in dim2
% w_corners_dim2 = center_dim2 ini_dim2 end_dim2
aux = 1:sq_size_pix:size_IMG_SEQ(2);
offset_center_dim2 = (size_IMG_SEQ(2)-aux(end))/2;
center_dim2 = offset_center_dim2 -(sq_dist_pix/2):sq_dist_pix:size_IMG_SEQ(2); 
ini_dim2 = round(center_dim2 - (sq_size_pix/2));
end_dim2 = round(center_dim2 + (sq_size_pix/2)-1);
center_dim2 = round(center_dim2);


% delete values out of the IMG_SEQ bounds (dim2)
ind2del = ini_dim2<1 | end_dim2>size_IMG_SEQ(2);
center_dim2(ind2del) = [];
ini_dim2(ind2del) = [];
end_dim2(ind2del) = [];


% generate all the posible center combinations
[center_dim1_2D,center_dim2_2D] = ndgrid(center_dim1,center_dim2);
[ini_dim1_2D,ini_dim2_2D] = ndgrid(ini_dim1,ini_dim2);
[end_dim1_2D,end_dim2_2D] = ndgrid(end_dim1,end_dim2);
N_fit_windows = numel(center_dim1_2D);

% resahape data to store coners in the following structure
% w_corners_dim1(1,:) = > centers dim1
% w_corners_dim1(2,:) = > initial corner dim1
% w_corners_dim1(3,:) = > end corner dim1
w_corners_dim1 = cat(1,reshape(center_dim1_2D,1,N_fit_windows),reshape(ini_dim1_2D,1,N_fit_windows),reshape(end_dim1_2D,1,N_fit_windows));

% w_corners_dim2(1,:) = > centers dim2
% w_corners_dim2(2,:) = > initial corner dim2
% w_corners_dim2(3,:) = > end corner dim2
w_corners_dim2 = cat(1,reshape(center_dim2_2D,1,N_fit_windows),reshape(ini_dim2_2D,1,N_fit_windows),reshape(end_dim2_2D,1,N_fit_windows));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% remove windows acccording to nan_percentage_thr  & mask_2D      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NAN_2D_percentage = nanmean(isnan(IMG_SEQ.IMG),3);
NAN_2D_mask = NAN_2D_percentage < nan_percentage_thr/100;
usable_data_mask_2D = mask_2D & NAN_2D_mask;

usable_data_percentage_win = nan(1,N_fit_windows);
average_depth_win = nan(1,N_fit_windows);

% check window by window the data percentage available & average depth
for i1 = 1:N_fit_windows

    % retrieve data according window limits
    data_usable_win = usable_data_mask_2D(w_corners_dim1(2,i1):w_corners_dim1(3,i1),w_corners_dim2(2,i1):w_corners_dim2(3,i1));

    % calculate percentage of usable data
    usable_data_percentage_win(i1) = nanmean(double(data_usable_win(:)));
    
    % retrieve the average depth for every window water_depth_mask_2D
    depth_win = water_depth_mask_2D(w_corners_dim1(2,i1):w_corners_dim1(3,i1),w_corners_dim2(2,i1):w_corners_dim2(3,i1));
    average_depth_win(i1) = nanmean(depth_win(:)); 
end

% delete windows with data under 100-nan_percentage_thr
ind2del = 100*usable_data_percentage_win < 100-nan_percentage_thr;
w_corners_dim1(:,ind2del) = [];
w_corners_dim2(:,ind2del) = [];
average_depth_win(ind2del) = [];
N_fit_windows = size(w_corners_dim1,2);
fitted_window_flag = false(1,N_fit_windows);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the General_data struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
image = IMG_SEQ.IMG(:,:,1);

Generic = struct('gridX',IMG_SEQ.gridX, ...
                'gridY',IMG_SEQ.gridY,'image',image,'Longitude',IMG_SEQ.Longitude,...
                'Latitude',IMG_SEQ.Latitude,'heading',IMG_SEQ.heading,...
                'altitude',IMG_SEQ.altitude,'time_stamp',IMG_SEQ.ts_video,...
                'usable_data_mask_2D',usable_data_mask_2D,...
                'water_depth_mask_2D',water_depth_mask_2D);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the Window struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Windows = struct('N_fit_windows',N_fit_windows,'w_corners_dim1',w_corners_dim1, ...
                'w_corners_dim2',w_corners_dim2,'sq_size_m', sq_size_m,...
                'sq_dist_m',sq_dist_m,'average_depth_win',average_depth_win,...
                'fitted_window_flag',fitted_window_flag);
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define K_steps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
if exist('K_steps','var') ==0 || isempty(K_steps)
    
    Nx = round(sq_size_m/IMG_SEQ.dx);
    Kx = (2 * pi * 1/IMG_SEQ.dx/Nx) * [-ceil((Nx-1)/2): floor((Nx-1)/2)];
    dkx = Kx(2)-Kx(1);

    pixel_size = sqrt(IMG_SEQ.dx.^2 + IMG_SEQ.dy.^2);
    dk_end = 2*pi/(4*pixel_size);
    dk_ini = 2*pi/(sq_size_m/4);
    
    step_ini = dk_ini:dkx:dk_end;
    step_end = step_ini + dkx;
    K_steps = cat(1,step_ini,step_end);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate current fit param struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create firist guess matrix
Ux_FG_vec =  Ux_limits_FG(1):U_FG_res:Ux_limits_FG(2);
Uy_FG_vec =  Uy_limits_FG(1):U_FG_res:Uy_limits_FG(2);
[Ux_FG_2D,Uy_FG_2D] = ndgrid(Ux_FG_vec,Uy_FG_vec);

% create second guess matrix (according max_current_shear_expected + U_SG_res)
Ux_SG_vec =  -max_current_shear_expected:U_SG_res:max_current_shear_expected;
Uy_SG_vec =  -max_current_shear_expected:U_SG_res:max_current_shear_expected;
[Ux_SG_2D,Uy_SG_2D] = ndgrid(Ux_SG_vec,Uy_SG_vec);

% get max min wavelengh
% waveLength_limits_m = 2*pi ./ [max(K_steps(:)) min(K_steps(:))];
K_limits = [2*pi/waveLength_limits_m(2) 2*pi/waveLength_limits_m(1)];
W_limits = [2*pi/wavePeriod_limits_sec(2) 2*pi/wavePeriod_limits_sec(1)];

% apply waveLength_limits_m to K_steps
% K_limits = [2*pi/waveLength_limits_m(2) 2*pi/waveLength_limits_m(1)];
ind2del_K_steps =  any(K_steps<K_limits(1) | K_steps>K_limits(2),1);
K_steps(:,ind2del_K_steps) = [];

fit_param = struct('Ux_FG_2D',Ux_FG_2D,'Uy_FG_2D',Uy_FG_2D,...
                   'Ux_SG_2D',Ux_SG_2D,'Uy_SG_2D',Uy_SG_2D,...
                   'w_width_FG',w_width_FG,'w_width_SG',w_width_SG,...
                   'waveLength_limits_m',waveLength_limits_m,...
                   'wavePeriod_limits_sec', wavePeriod_limits_sec,...
                   'K_limits',K_limits,'W_limits',W_limits,...
                   'K_steps',K_steps,...
                   'max_current_shear_expected',max_current_shear_expected);
               

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate current fit output struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%               

% TO D0
out_fit = [];


% %%%%%%%%%%%%%%%%%%
% % generate STCFIT
% %%%%%%%%%%%%%%%%%
STCFIT_multi_K = struct('Generic',Generic,'Windows',Windows,'fit_param',fit_param,'out_fit',out_fit);




end
