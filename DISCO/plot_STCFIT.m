function [h] = plot_STCFIT(STCFIT)
% plot STCFIT. Usefull to check the windows distribution.
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


% plot image
h = figure('units','normalized','outerposition',[0 0 1 1]);

pcolor(STCFIT.Generic.gridX, STCFIT.Generic.gridY, STCFIT.Generic.image);
shading flat;
axis xy equal tight;
hold on;
colormap(gray);

% plot window centers
IND_center = sub2ind(size(STCFIT.Generic.gridX),STCFIT.Windows.w_corners_dim1(1,:),STCFIT.Windows.w_corners_dim2(1,:));
centers_x = STCFIT.Generic.gridX(IND_center);
centers_y = STCFIT.Generic.gridY(IND_center);
plot(centers_x,centers_y,'+r','Linewidth',1.5);

offset_x_text = STCFIT.Windows.sq_size_m/8;

% plot window
for i1 = 1:STCFIT.Windows.N_fit_windows

    % plot square
    x_ini =STCFIT.Generic.gridX(STCFIT.Windows.w_corners_dim1(2,i1),STCFIT.Windows.w_corners_dim2(2,i1));
    x_end =STCFIT.Generic.gridX(STCFIT.Windows.w_corners_dim1(3,i1),STCFIT.Windows.w_corners_dim2(3,i1));
    y_ini =STCFIT.Generic.gridY(STCFIT.Windows.w_corners_dim1(2,i1),STCFIT.Windows.w_corners_dim2(2,i1));
    y_end =STCFIT.Generic.gridY(STCFIT.Windows.w_corners_dim1(3,i1),STCFIT.Windows.w_corners_dim2(3,i1));
    
    color_mat = rand(1,3);
    plot([x_ini x_ini x_end x_end x_ini], [y_ini y_end y_end y_ini y_ini],'Color',color_mat,'LineStyle','--');

    % plot square number
    x_text = STCFIT.Generic.gridX(STCFIT.Windows.w_corners_dim1(1,i1),STCFIT.Windows.w_corners_dim2(1,i1)) + offset_x_text;
    y_text = STCFIT.Generic.gridY(STCFIT.Windows.w_corners_dim1(1,i1),STCFIT.Windows.w_corners_dim2(1,i1));
    text(x_text,y_text,num2str(i1),'HorizontalAlignment','left','Color',color_mat);
end

% % plot output if is avaible
% if ~isempty(STCFIT.out_fit)
%     out_fit = [STCFIT.out_fit(:).SG_fit];
%     Ux_fit = [out_fit(:).Ux_fit];
%     Uy_fit = [out_fit(:).Uy_fit];
%     SNR_max = [out_fit(:).SNR_max];
%     SNR_density_max = [out_fit(:).SNR_density_max];
%     
%     % filter with SNR
%     ind2use = SNR_density_max>3;
%     
%     % plot 
%     quiver(centers_x(ind2use),centers_y(ind2use),Ux_fit(ind2use),Uy_fit(ind2use),'k');
%     % U_mag = sqrt(Ux_fit(ind2use).^2 + Uy_fit(ind2use).^2);
%     
% end

end

