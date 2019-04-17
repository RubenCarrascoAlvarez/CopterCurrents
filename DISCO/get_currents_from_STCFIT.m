function [UTM_currents, Camera_currents] = get_currents_from_STCFIT(STCFIT,SNR_thr,SNR_density_thr,currentdir_flag,reprocess_SNR_flag)
% Get current from STCFIT structure. Very usefull for plotting currents.
%  Note 1: 'run_current_fit.m' must be run before in the STCFIT structure.
%  Note 2: UTM_currents requires the deg2utm external libray.
%
% Input:
%   STCFIT: fit structure
%   SNR_thr: signal to noise ratio threshold for valid data. [0.2]
%   SNR_density_thr:  signal to noise density ratio threshold for valid data [5].
%   currentdir_flag: [1] current direction 'going to'
%                     0  current direction 'comming from'
%   reprocess_SNR_flag: processs SNR with SNR_thr and SNR_density_thr
%                       thresholds. Only the second guess is processed. 
%                        1   reprocess.
%                       [0]  not reprocess.
%
% Output:
%   UTM_currents Structure: (The reference is UTM zone grid)
%     UTM_currents.x: vector with UTM x postion
%     UTM_currents.y: vector with UTM y postion
% 	  UTM_currents.Ux: vector with Ux in meter/seconds
% 	  UTM_currents.Uy: vector with Uy in meter/seconds
% 	  UTM_currents.SNR: SNR vector
%     UTM_currents.SNR_density: SNR density vector
% 	  UTM_currents.time_stamp: time stamp in datenum format
% 	  UTM_currents.utmzone: ut'32 U'
%     UTM_currents.uav_x_utm: x drone position in utm coordenates.
%     UTM_currents.uav_y_utm: y drone position in utm coordenates.
%     UTM_currents.currentdir_flag: [1] current direction 'going to'
%                                   [0] current direction 'comming from'
%
%   Camera_currents Structure: (The reference is in the center of the camera)
%     Camera_currents.x: x postion from camera center
%     Camera_currents.y: y postion from camera center
% 	  Camera_currents.Ux:  Ux [m/s] in the camera direction.
% 	  Camera_currents.Uy:  Uy [m/s] in the camera direction.
% 	  Camera_currents.SNR: SNR vector
%     Camera_currents.SNR_density: SNR density vector
% 	  Camera_currents.time_stamp: time stamp in datenum format
%     Camera_currents.currentdir_flag: [1] current direction 'going to'
%                                      [0] current direction 'comming from'
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% (C) 2019, Ruben Carrasco <ruben.carrasco@hzg.de>
%
% This file is part of DISCO
%
% DISCO has been developed by Department of Radar Hydrography, at
% Helmholtz-Zentrum Geesthacht Centre for Materials and Coastal Research 
% (Germany), based in the work exposed in: 
% 
% M. Streßer, R. Carrasco and J. Horstmann, "Video-Based Estimation 
% of Surface Currents Using a Low-Cost Quadcopter," in IEEE Geoscience 
% and Remote Sensing Letters, vol. 14, no. 11, pp. 2027-2031, Nov. 2017.
% doi: 10.1109/LGRS.2017.2749120
%
%
% DISCO is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% DISCO is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with DISCO.  If not, see <http://www.gnu.org/licenses/>.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if exist('SNR_thr','var') == 0 || isempty(SNR_thr)
    SNR_thr = 0.2;
end

if exist('SNR_density_thr','var') == 0 || isempty(SNR_density_thr)
    SNR_density_thr = 5;
end

if exist('currentdir_flag','var') == 0 || isempty(currentdir_flag)
    currentdir_flag = 1;
end

if exist('reprocess_SNR_flag','var') == 0 || isempty(reprocess_SNR_flag)
    reprocess_SNR_flag = 0;
end




% check output
if isempty(STCFIT.out_fit)
   error('get_currents_from_STCFIT: the structure STCFIT was not processed. Run run_current_fit.m before.'); 
end

if reprocess_SNR_flag == 1
    % reprocess fit
    for i1 = 1:length(STCFIT.out_fit)
         % [Ux,Uy,SNR_max,lin_ind] = retrieve_best_fit_from_SNR_2D(SNR_2D,Ux_2D,Uy_2D,method,display)
         [STCFIT.out_fit(i1).SG_fit.Ux_fit,...
             STCFIT.out_fit(i1).SG_fit.Uy_fit,...
             STCFIT.out_fit(i1).SG_fit.SNR_density_max,...
             lin_ind] = retrieve_best_fit_from_SNR_2D(...
             STCFIT.out_fit(i1).SG_fit.SNR_density_2D,...
             STCFIT.out_fit(i1).SG_fit.Ux_2D,...
             STCFIT.out_fit(i1).SG_fit.Uy_2D,0,0);
            
          STCFIT.out_fit(i1).SG_fit.SNR_max =  STCFIT.out_fit(i1).SG_fit.SNR_2D(lin_ind);
        
    end
    
end


% get signal to noise data 
SG_fit  = [STCFIT.out_fit(:).SG_fit]; 
SNR = [SG_fit(:).SNR_max];
SNR_density = [SG_fit(:).SNR_density_max];

% get index to set to NaN (low SNR cases)
ind2setnan = (SNR < SNR_thr) | (SNR_density  < SNR_density_thr);

% retrieve currents from second fit
Ux = [SG_fit(:).Ux_fit];
Ux(ind2setnan) = NaN;
Uy = [SG_fit(:).Uy_fit];
Uy(ind2setnan) = NaN;

% exchange direction
if currentdir_flag == 1
    Ux = -Ux;
    Uy = -Uy;
end

% get position from camera
IND_center = sub2ind(size(STCFIT.Generic.gridX),STCFIT.Windows.w_corners_dim1(1,:),STCFIT.Windows.w_corners_dim2(1,:));
x = STCFIT.Generic.gridX(IND_center);
y = STCFIT.Generic.gridY(IND_center);

% generate camera strucure
Camera_currents = struct('x',x,'y',y,'Ux',Ux,'Uy',Uy,'SNR',SNR,...
    'SNR_density',SNR_density,'time_stamp',STCFIT.Generic.time_stamp,...
    'currentdir_flag',currentdir_flag);
              

% convert to UTM      
try
    [uav_x_utm,uav_y_utm,utmzone] = deg2utm(STCFIT.Generic.Latitude,STCFIT.Generic.Longitude);
    % rotate x and y coordenates according to heading
    [grid_center_pts] = rotatePoint3d([x(:) y(:) zeros(length(x(:)),1)],STCFIT.Generic.heading,0,0); % yaw from video exif
    x_utm = reshape(grid_center_pts(:,1),size(x)) + uav_x_utm;
    y_utm = reshape(grid_center_pts(:,2),size(x)) + uav_y_utm;   

    % rotate Ux Uy
    [grid_U_rot] = rotatePoint3d([Ux(:) Uy(:) zeros(length(Ux(:)),1)],STCFIT.Generic.heading,0,0);
    Ux_utm = reshape(grid_U_rot(:,1),size(Ux));
    Uy_utm = reshape(grid_U_rot(:,2),size(Uy));
    
    % generate UTM strucure
    UTM_currents = struct('x',x_utm,'y',y_utm,'Ux',Ux_utm,'Uy',Uy_utm,'SNR',SNR,...
        'SNR_density',SNR_density,'time_stamp',STCFIT.Generic.time_stamp,...
        'utmzone',utmzone,'uav_x_utm',uav_x_utm,'uav_y_utm',uav_y_utm, ...
        'currentdir_flag',currentdir_flag);

    
catch
    disp('Check if the external library deg2utm was installed');
    warning('get_currents_from_STCFIT: UTM coordenates not retrieved.');
    
    UTM_currents = [];
end



end

