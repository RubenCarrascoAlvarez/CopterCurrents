function folder = generate_SNR_plots_from_SCTFIT(STCFIT,folder)
% Generate pictures with all SNR plots from STCFIT
%
%   Input:
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
%
%         STCFIT.Windows.N_fit_windows: number of windows used for the current fit.
%         STCFIT.Windows.w_corners_dim1: 3xN_fit_windows
%         STCFIT.Windows.w_corners_dim2: 3xN_fit_windows
%         STCFIT.Windows.sq_size_m: square size in meter for the current fit windows.
%         STCFIT.Windows.sq_dist_m: square distance in meter between current fit windows.
%
%       folder: folder to save the plots
%
%   Output:
%       folder: folder to save the plots
%
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




if exist('folder','var') == 0
    
    if isunix 
        folder = [pwd '/' 'output/']; 
    else
        folder = [pwd '\' 'output\'];
    end
end

if exist(folder,'dir') == 0
    mkdir(folder); 
    disp('Output folder created.')
end


for n_windowfitted = 1:STCFIT.Windows.N_fit_windows
   
    % generate plot
    [h] = plot_window_fit_result(STCFIT,n_windowfitted);
    
    % save plot
    savepng = [folder 'SNR_' num2str(n_windowfitted,'%04.0f') '.png'];
    print(h,'-dpng',savepng);
    
    % close plot
    close(h);
    
end


end