function [h] = plot_window_fit_result(STCFIT,n_windowfitted)
% plot Signal to Noise maps result on the current fit.  
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
%   n_windowfitted: number of window to plot
%
%   Output:
%       h =  handle to figure
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

Linewidth = 3;

if n_windowfitted>STCFIT.Windows.N_fit_windows
   error(['plot_window_fit_result: n_windowfitted ' num2str(n_windowfitted) ' does not valid']) 
end

if isempty(STCFIT.out_fit)
   error('plot_window_fit_result: STCFIT struct not fitted yet')  
   
else
   out_fit = STCFIT.out_fit(n_windowfitted);
   
   % create figure
   name_str = ['Window fit ' num2str(STCFIT.Windows.N_fit_windows)];
   h = figure('units','normalized','outerposition',[0 0 1 1],'name',name_str);
   colormap(jet);
   
   offset_x =  (out_fit.FG_fit.Ux_2D(2,1) - out_fit.FG_fit.Ux_2D(1,1))/2;
   offset_y =  (out_fit.FG_fit.Uy_2D(1,2) - out_fit.FG_fit.Uy_2D(1,1))/2;
   
   % plot first guess
   str_U_FG = ['Ux = ' num2str(out_fit.FG_fit.Ux_fit) ' [m/s]'  '   Uy = ' num2str(out_fit.FG_fit.Uy_fit) ' [m/s]'];

   ax(1) = subplot(2,2,1);
   pcolor(out_fit.FG_fit.Ux_2D - offset_x, out_fit.FG_fit.Uy_2D - offset_y, out_fit.FG_fit.SNR_2D);
   shading flat;
   axis xy equal tight;
   hold on;
   plot(out_fit.FG_fit.Ux_fit, out_fit.FG_fit.Uy_fit,'+k','Linewidth',Linewidth);
   xlabel('Ux [m/s]');
   ylabel('Uy [m/s]');
   cb1 = colorbar;
   ylabel(cb1,'SNR')
   title({['First guess: ' str_U_FG ],['SNR = ' num2str(out_fit.FG_fit.SNR_max)]});

   
   ax(2) = subplot(2,2,2);
   pcolor(out_fit.FG_fit.Ux_2D - offset_x, out_fit.FG_fit.Uy_2D - offset_y, out_fit.FG_fit.SNR_density_2D);
   shading flat;
   axis xy equal tight;
   hold on;
   plot(out_fit.FG_fit.Ux_fit, out_fit.FG_fit.Uy_fit,'+k','Linewidth',Linewidth);
   xlabel('Ux [m/s]');
   ylabel('Uy [m/s]');
   cb2 = colorbar;
   ylabel(cb2,'SNR density')
   title({['First guess: ' str_U_FG ],['SNR density = ' num2str(out_fit.FG_fit.SNR_density_max)]});
   
   
   linkaxes(ax(1:2),'xy');
   
   % plot second guess 
   
   offset_x =  (out_fit.SG_fit.Ux_2D(2,1) - out_fit.SG_fit.Ux_2D(1,1))/2;
   offset_y =  (out_fit.SG_fit.Uy_2D(1,2) - out_fit.SG_fit.Uy_2D(1,1))/2;
   
   str_U_SG = ['Ux = ' num2str(out_fit.SG_fit.Ux_fit) ' [m/s]'  '   Uy = ' num2str(out_fit.SG_fit.Uy_fit) ' [m/s]'];
   ax(3) = subplot(2,2,3);
   pcolor(out_fit.SG_fit.Ux_2D - offset_x, out_fit.SG_fit.Uy_2D - offset_y, out_fit.SG_fit.SNR_2D);
   shading flat;
   axis xy equal tight;
   hold on;
   plot(out_fit.SG_fit.Ux_fit, out_fit.SG_fit.Uy_fit,'+k','Linewidth',Linewidth);
   xlabel('Ux [m/s]');
   ylabel('Uy [m/s]');
   cb3 = colorbar;
   ylabel(cb3,'SNR')
   title({['Second guess: ' str_U_SG ],['SNR = ' num2str(out_fit.SG_fit.SNR_max)]});

   
   ax(4) = subplot(2,2,4);
   pcolor(out_fit.SG_fit.Ux_2D - offset_x, out_fit.SG_fit.Uy_2D - offset_y, out_fit.SG_fit.SNR_density_2D);
   shading flat;
   axis xy equal tight;
   hold on;
   plot(out_fit.SG_fit.Ux_fit, out_fit.SG_fit.Uy_fit,'+k','Linewidth',Linewidth);
   xlabel('Ux [m/s]');
   ylabel('Uy [m/s]');
   cb4 = colorbar;
   ylabel(cb4,'SNR density')
   title({['Second guess: ' str_U_SG ],['SNR density = ' num2str(out_fit.SG_fit.SNR_density_max)]});
   
   linkaxes(ax(3:4),'xy');
   
end



end

