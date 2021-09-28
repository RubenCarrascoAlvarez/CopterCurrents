function STCFIT = run_current_fit_depth_profile(IMG_SEQ,STCFIT,windowList)

% Run wave dispersion realtion current fit, in every window defined in the STCFIT
%   
%   Input:
%     IMG_SEQ structure
%         IMG_SEQ.IMG: 3D Intensty data (Nx,Ny,Nt)
%         IMG_SEQ.gridX: 2D grid X data in meters (Nx,Ny)
%         IMG_SEQ.gridY: 2D grid Y data in meters (Nx,Ny)
%         IMG_SEQ.dx: X resolution in meters.
%         IMG_SEQ.dy: Y resolution in meters.
%         IMG_SEQ.dt: time resolution in seconds.
%         IMG_SEQ.ts_video: video time stamps in datenum format. 1xNt vector.
%         IMG_SEQ.altitude: altitude to the water surface in meters.
%         IMG_SEQ.pitch: video pitch (for Nadir pitch = -90)
%         IMG_SEQ.roll: video roll (for Nadir roll = 0)
%         IMG_SEQ.heading: video heading to North (yaw)
%         IMG_SEQ.Longitude: video longitude in degrees
%         IMG_SEQ.Latitude: video latitude in degrees
%         IMG_SEQ.mediainfo_string: raw mediainfo string
%         IMG_SEQ.Georeference_Struct_config: Georeference_Struct used
%
%     STCFIT: Structure defining the current fit input parameters
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
%         STCFIT.Windows.average_depth_win: average depth in meters.
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
%   STCFIT.out_fit with the following fields: Where FG_fit (First Guess) contains:
%      FG_fit.Ux_2D: 2d matrix with x first guess 
%      FG_fit.Uy_2D: 2d matrix with y first guess 
%      FG_fit.signal_2D: Signal in the dispersion relation in every (Ux,Uy) pair
%      FG_fit.noise_2D: Noise out of the dispersion relation in every (Ux,Uy) pair
%      FG_fit.signal_nvalues_2D: n values in the dispersion relation in every (Ux,Uy) pair
%      FG_fit.noise_nvalues_2D: n values out of the dispersion relation in every (Ux,Uy) pair
%      FG_fit.SNR_2D: Signal to Noise ratio in every (Ux,Uy) pair
%      FG_fit.SNR_density_2D: Signal to Noise density ratio in every (Ux,Uy) pair
%      FG_fit.Ux_fit: best first guess fit for Ux
%      FG_fit.Uy_fit: best first guess fit for Uy
%      FG_fit.SNR_max: SNR corresponding to the best first guess
%      FG_fit.SNR_density_max: SNR density corresponding to the best first guess
%    
%   Where SG_fit (Second Guess) contains:
%      SG_fit.Ux_2D: 2d matrix with x second guess 
%      SG_fit.Uy_2D: 2d matrix with y second guess 
%      SG_fit.signal_2D: Signal in the dispersion relation in every (Ux,Uy) pair
%      SG_fit.noise_2D: Noise out of the dispersion relation in every (Ux,Uy) pair
%      SG_fit.signal_nvalues_2D: n values in the dispersion relation in every (Ux,Uy) pair
%      SG_fit.noise_nvalues_2D: n values out of the dispersion relation in every (Ux,Uy) pair
%      SG_fit.SNR_2D: Signal to Noise ratio in every (Ux,Uy) pair
%      SG_fit.SNR_density_2D: Signal to Noise density ratio in every (Ux,Uy) pair
%      SG_fit.Ux_fit: best second guess fit for Ux
%      SG_fit.Uy_fit: best second guess fit for Uy
%      SG_fit.SNR_max: SNR corresponding to the best second guess
%      SG_fit.SNR_density_max: SNR density corresponding to the second first guess
%    
% Note the final result of every window is stored in
%      STCFIT.out_fit(N)SG_fit.Ux_fit: X current in meters per second
%      STCFIT.out_fit(N)SG_fit.Uy_fit: Y current in meters per second
%       where N is the window number and goes from 1 to STCFIT.Windows.N_fit_windows
%

%Additonal STCFIT.NSPP_fit_param: 
%       
%       fit_param.wavenumbers: list of wavenumbers in [rad/length] at which Doppler shift velocities will be extracted
%       fit_param.include2ndHarmonic (optional): whether to include 2nd harmonic of the spectrum in the fit (set to false by default)
%       fit_param.logFlag (optional): whether to do the fit in log space (false by default)
%       fit_param.omegaWidthFun: function handle as a function of wavenumber i.e.
%           @(k) f(k)...., specifying frequency width of the weighting function in
%           frequency-angle space (constant wavenumber). Width is half-width 1/e^2
%           point of a Gaussian function

%The following parameters involve post-processing of the Doppler shifts:
%       fit_param.SNR_filter: whether to use a signal-to-noise filter (false by
%       fit_param.SNR_threshold: threshold signal-to-noise value for above filter (set to 2.0 by default)
%       fit_param.Peak_filter: whether to use a multiple peaks filter (false by default)
%       fit_param.Peak_threshold: peak threshold of maximum value (0.5 by default)
%       fit_param.Outlier_filter: whether to use an outlier filter (quartile-based) (false by default)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% (C) 2019, Ruben Carrasco <ruben.carrasco@hzg.de> and (C) 2020, Benjamin
% Smeltzer (benjamin.smeltzer@ntnu.no)
%
% This file is part of CopterCurrents
%
% CopterCurrents has been developed by Department of Radar Hydrography, at
% Helmholtz-Zentrum Geesthacht Centre for Materials and Coastal Research 
% (Germany), based in the work exposed in: 
% 
% M. StreÃŸer, R. Carrasco and J. Horstmann, "Video-Based Estimation 
% of Surface Currents Using a Low-Cost Quadcopter," in IEEE Geoscience 
% and Remote Sensing Letters, vol. 14, no. 11, pp. 2027-2031, Nov. 2017.
% doi: 10.1109/LGRS.2017.2749120

%The PEDM section is described in the following article:

%Smeltzer, B. K., Æsøy, E., Ådnøy, A.,& Ellingsen, S. Å. (2019). An improved
%method for determining near-surface currents from wave dispersion measurements. 
%Journal of Geophysical Research: Oceans, 124. https://doi.org/10.1029/2019JC015202


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

if ~exist('windowList','var')
    windowList = 1:STCFIT.Windows.N_fit_windows;
end

if isempty(windowList)
    windowList = 1:STCFIT.Windows.N_fit_windows;
end
    

g = 9.81;%Gravitatinal acceleration
T = 0.072 / 1000;%Surace tension / density coefficient.

for iWin = 1:numel(windowList)
    i1 = windowList(iWin);
    
    % display window fit number
    disp(['Window ' num2str(i1) ' of ' num2str(STCFIT.Windows.N_fit_windows)]);
    
    % cut i1 window
    IMG_SEQ_Window = double(IMG_SEQ.IMG(STCFIT.Windows.w_corners_dim1(2,i1):STCFIT.Windows.w_corners_dim1(3,i1),...
                                        STCFIT.Windows.w_corners_dim2(2,i1):STCFIT.Windows.w_corners_dim2(3,i1),:));
                                    
    % get spectrum
    Spectrum = retrieve_power_spectrum(IMG_SEQ_Window,IMG_SEQ.dx, IMG_SEQ.dy, IMG_SEQ.dt,STCFIT.fit_param.K_limits, STCFIT.fit_param.W_limits); 
    
   
    fit_param = STCFIT.NSPP_fit_param;

   
    fit_param.Ux_2D = STCFIT.fit_param.Ux_FG_2D;
    fit_param.Uy_2D = STCFIT.fit_param.Uy_FG_2D;

    
    h = STCFIT.Windows.average_depth_win(i1);%Local water depth.
    
    Properties = struct('h',h,'g',g,'T',T);
    
    %Get Doppler shift velocities
    %tic;
    out_DS_i = get_doppler_shift_velocities_nsp(Spectrum,fit_param,Properties);
    %toc;


out(iWin) = find_current_depth_profile(out_DS_i.wavenumbers,out_DS_i.Ux_filt,out_DS_i.Uy_filt,h);
out(iWin).NSP = out_DS_i;
end


STCFIT.out_depth_profile = out;


end

