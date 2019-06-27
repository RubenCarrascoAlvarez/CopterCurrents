function h  = plot_currents_map(currents,STCFIT,arrow_scale)
% plot current map
%
% Inputs:
%   currents: current structure from get_currents_from_STCFIT. 
%       UTM_currents or Camera_currents current structure.
%   STCFIT: STCFIT structure
%   arrow_scale: plotting arrows scale in meters.
%
% Outputs:
%   h = handle to the figure
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

% get grid
gridX = STCFIT.Generic.gridX;
gridY = STCFIT.Generic.gridY;

if isfield(currents,'utmzone') == 1 % UTM current stucture
    
    % rotate x and y coordenates according to heading

    [grid_center_pts] = rotatePoint3d([gridX(:) gridY(:) zeros(length(gridX(:)),1)],STCFIT.Generic.heading,0,0); 
    gridX = reshape(grid_center_pts(:,1),size(gridX)) + currents.uav_x_utm;
    gridY = reshape(grid_center_pts(:,2),size(gridX)) + currents.uav_y_utm; 
    
    title_str = ['UTM zone: ' currents.utmzone '  date: ' datestr(currents.time_stamp(1))];
else
    title_str = ['Orientation to North: ' num2str(STCFIT.Generic.heading) ' [deg].  date: ' datestr(currents.time_stamp(1))];
    
end

% get scale factor for quiver
if exist('arrow_scale','var') == 0
    grid_diag_size = sqrt((max(gridX(:)) - min(gridX(:)))^2  + (max(gridY(:)) - min(gridY(:)))^2);
    arrow_scale = grid_diag_size/10;
end


% plot
h = figure('units','normalized','outerposition',[0 0 1 1]);
colormap(gray);
pcolor(gridX, gridY, STCFIT.Generic.image);
shading flat;
axis xy equal tight;
hold on;
quiver(currents.x,currents.y,currents.Ux *arrow_scale,currents.Uy*arrow_scale,...
       'b','Linewidth',3,'AutoScale','off');
title(title_str);   
   
if isfield(currents,'utmzone') == 1 % UTM current stucture
    xlabel(['UTM X grid ' currents.utmzone ' [meters]']);
    ylabel(['UTM Y grid ' currents.utmzone ' [meters]']);
else
    xlabel('X Distance to camera center [meters]');
    ylabel('Y Distance to camera center [meters]');
end


end

