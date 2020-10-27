function fit_out = fit_Spectrum2dispersionRelation_multi_K(Spectrum,fit_param,water_depth, window_used)
% fit Spectrum to the waves dispersion relation fundamental mode.
%    The fit is made in 2 steps, first guess and second. Both step data are
%    store in fit_out.
%
%    Spectrum.power_Spectrum: (Kx,Ky,W) power Spectrum
%    Spectrum.Kx_3D: 3D Kx grid corresponding to Spectrum.power_Spectrum [rad/m]
%    Spectrum.Ky_3D: 3D Ky grid corresponding to Spectrum.power_Spectrum [rad/m]
%    Spectrum.W_3D: 3D W grid corresponding to Spectrum.power_Spectrum [rad/sec]
%    Spectrum.dKx: Kx resolution [rad/m]
%    Spectrum.dKy: Ky resolution [rad/m]
%    Spectrum.dW: W resolution [rad/sec]
%    Spectrum.Kx_orig_limits: raw spectra Kx limits (before the cut) [rad/m]
%    Spectrum.Ky_orig_limits: raw spectra Ky limits (before the cut) [rad/m]
%    Spectrum.W_orig_limits: raw spectra W limits (before the cut) [rad/sec]
%
%    STCFIT.fit_param structure
%    fit_param.Ux_FG_2D: MxN Ux first guess matrix
%    fit_param.Uy_FG_2D: MxN Uy first guess matrix
%    fit_param.Ux_SG_2D: MxN Ux second guess matrix (offset matrix)
%    fit_param.Uy_SG_2D: MxN Uy second guess matrix (offset matrix) 
%    fit_param.w_width_FG: first guess filter width in w [rad/s] 
%    fit_param.w_width_SG: second guess filter width in w [rad/s]
%    fit_param.waveLength_limits_m: [min max] wavelength to use [meters]
%    fit_param.wavePeriod_limits_sec: [min max] wave Period to use [seconds]
%    fit_param.K_limits: [min max] wave number to use [rad/m]
%    fit_param.W_limits: [min max] wave frequency to use [rad/sec]
%
%    water_depth: water depth in meters.
%    window_used: number of the window used.
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  First Guess
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create vectors
Ux_FG_vec = fit_param.Ux_FG_2D(:);
Uy_FG_vec = fit_param.Uy_FG_2D(:);
w_width =   fit_param.w_width_FG;

signal_vec = nan(size(Ux_FG_vec));
signal_nvalues = nan(size(Ux_FG_vec));

noise_vec = nan(size(Ux_FG_vec));
noise_nvalues = nan(size(Ux_FG_vec));

% for i1 = 1:length(Ux_FG_vec)
parfor i1 = 1:length(Ux_FG_vec)

    % fit dispersion relation
    [DS_3D_mask] = get_Dispersion_Relation_3D_mask(Spectrum, water_depth, Ux_FG_vec(i1), Uy_FG_vec(i1), w_width);
    
    % get signal to noise ratio
    signal = Spectrum.power_Spectrum(DS_3D_mask);
    noise = Spectrum.power_Spectrum(DS_3D_mask==0);
    
    signal_vec(i1) = nansum(signal);
    noise_vec(i1) = nansum(noise);
    
    signal_nvalues(i1) = sum(isfinite(signal));
    noise_nvalues(i1) = sum(isfinite(noise));
    
end

% resahpe out vectors
signal_FG_2D = reshape(signal_vec,size(fit_param.Ux_FG_2D));
noise_FG_2D = reshape(noise_vec,size(fit_param.Ux_FG_2D));

signal_nvalues_2D = reshape(signal_nvalues,size(fit_param.Ux_FG_2D));
noise_nvalues_2D = reshape(noise_nvalues,size(fit_param.Ux_FG_2D));

% get signal to noise ratio
SNR_FG =  signal_FG_2D./noise_FG_2D;
SNR_density_FG =  (signal_FG_2D./ signal_nvalues_2D) ./ (noise_FG_2D ./ noise_nvalues_2D);

% retrieve best fit from SNR_density
[Ux_fit_FG,Uy_fit_FG,~,lin_ind] = retrieve_best_fit_from_SNR_2D(SNR_density_FG,fit_param.Ux_FG_2D,fit_param.Uy_FG_2D,0,0);
SNR_FG_max = SNR_FG(lin_ind);
SNR_density_FG_max = SNR_density_FG(lin_ind);

% display first guess result
disp(['Ux_Fg: ' num2str(Ux_fit_FG) ' m/s   Uy_Fg: ' num2str(Uy_fit_FG) ' m/s   SNR density: ' num2str(SNR_density_FG_max) ]);

% create first guess sctructure
FG_fit =  struct('signal_2D',signal_FG_2D,'noise_2D',noise_FG_2D,...
       'signal_nvalues_2D',signal_nvalues_2D,'noise_nvalues_2D',noise_nvalues_2D,...
       'SNR_2D',SNR_FG,'SNR_density_2D',SNR_density_FG,...
       'Ux_fit',Ux_fit_FG,'Uy_fit',Uy_fit_FG,...
       'SNR_max',SNR_FG_max,'SNR_density_max',SNR_density_FG_max,...
       'Ux_2D',fit_param.Ux_FG_2D,'Uy_2D',fit_param.Uy_FG_2D);

% note: the first guess is processed as before.
% the multi K fit is applied only for the second guess



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Second Guess
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% get K from the spectrum data
K = sqrt(Spectrum.Kx_3D.^2 + Spectrum.Ky_3D.^2);

% get index valid in the spectrum according to fit_param.K_steps
index_valid_K_steps = cell(1,size(fit_param.K_steps,2));

for i2 = 1:size(fit_param.K_steps,2)
    index_valid_K_steps{i2} =  K >= fit_param.K_steps(1,i2) & ...
                               K <= fit_param.K_steps(2,i2) ;                      
%     index_valid_K_steps{i2} =  find(K >= fit_param.K_steps(1,i2) & ...
%                                     K <= fit_param.K_steps(2,i2)    );                       
end


% create SG vectors
Ux_SG_2D = fit_param.Ux_SG_2D + Ux_fit_FG;
Uy_SG_2D = fit_param.Uy_SG_2D + Uy_fit_FG;

Ux_SG_vec = Ux_SG_2D(:);
Uy_SG_vec = Uy_SG_2D(:);
w_width =   fit_param.w_width_SG;

signal_vec = nan(size(Ux_SG_vec,1),size(fit_param.K_steps,2));
signal_nvalues = nan(size(Ux_SG_vec,1),size(fit_param.K_steps,2));

noise_vec = nan(size(Ux_SG_vec,1),size(fit_param.K_steps,2));
noise_nvalues = nan(size(Ux_SG_vec,1),size(fit_param.K_steps,2));

K_steps_length = size(fit_param.K_steps,2);



% for i1 = 1:length(Ux_SG_vec)
parfor i1 = 1:length(Ux_SG_vec)
    
    disp(['second gues ' num2str(i1) ' of ' num2str(length(Ux_SG_vec))]);
    
    % fit dispersion relation
    [DS_3D_mask] = get_Dispersion_Relation_3D_mask(Spectrum, water_depth, Ux_SG_vec(i1), Uy_SG_vec(i1), w_width);
    
    % get signal to noise ratio
    for i2 = 1:K_steps_length
        
        % cut DS_3D_mask according kstep
        DS_3D_mask_signal = DS_3D_mask & index_valid_K_steps{i2};
        DS_3D_mask_noise = ~DS_3D_mask & index_valid_K_steps{i2};
        
        signal = Spectrum.power_Spectrum(DS_3D_mask_signal);
        noise = Spectrum.power_Spectrum(DS_3D_mask_noise);

        signal_vec(i1,i2) = nansum(signal);
        noise_vec(i1,i2) = nansum(noise);

        signal_nvalues(i1,i2) = sum(isfinite(signal));
        noise_nvalues(i1,i2) = sum(isfinite(noise));
    end
    
end

% resahpe out vectors

% signal_SG_3D = nan(size(fit_param.Ux_SG_2D,1),size(fit_param.Ux_SG_2D,2),K_steps_length);
% noise_SG_3D = nan(size(fit_param.Ux_SG_2D,1),size(fit_param.Ux_SG_2D,2),K_steps_length);
% signal_nvalues_3D  = nan(size(fit_param.Ux_SG_2D,1),size(fit_param.Ux_SG_2D,2),K_steps_length);
% noise_nvalues_3D  = nan(size(fit_param.Ux_SG_2D,1),size(fit_param.Ux_SG_2D,2),K_steps_length);
% 
% for i2 = 1:K_steps_length
%     signal_SG_3D(:,:,i2) = reshape(signal_vec(:,i2),size(fit_param.Ux_SG_2D));
%     noise_SG_3D(:,:,i2) = reshape(noise_vec(:,i2),size(fit_param.Ux_SG_2D));
%     signal_nvalues_3D(:,:,i2) = reshape(signal_nvalues(:,i2),size(fit_param.Ux_SG_2D));
%     noise_nvalues_3D(:,:,i2) = reshape(noise_nvalues(:,i2),size(fit_param.Ux_SG_2D));
% end


signal_SG_3D = reshape(signal_vec,[size(fit_param.Ux_SG_2D,1) size(fit_param.Ux_SG_2D,2) K_steps_length]);
noise_SG_3D = reshape(noise_vec,[size(fit_param.Ux_SG_2D,1) size(fit_param.Ux_SG_2D,2) K_steps_length]);
signal_nvalues_3D = reshape(signal_nvalues,[size(fit_param.Ux_SG_2D,1) size(fit_param.Ux_SG_2D,2) K_steps_length]);
noise_nvalues_3D = reshape(noise_nvalues,[size(fit_param.Ux_SG_2D,1) size(fit_param.Ux_SG_2D,2) K_steps_length]);


% get signal to noise ratio
SNR_SG =  signal_SG_3D./noise_SG_3D;
SNR_density_SG =  (signal_SG_3D./ signal_nvalues_3D) ./ (noise_SG_3D ./ noise_nvalues_3D);


% retrieve best fit from SNR_density in every depth
Ux_fit_SG = nan(1,K_steps_length);
Uy_fit_SG = nan(1,K_steps_length);
SNR_SG_max = nan(1,K_steps_length);
SNR_density_SG_max = nan(1,K_steps_length);
% SNR_SG_dispersion = nan(1,K_steps_length);
% SNR_density_SG_dispersion = nan(1,K_steps_length);

K_steps_AV = mean(fit_param.K_steps,1);

% dispersion_conf_area_mps = 0.05;

for i2 = 1:K_steps_length
    
    SNR_density_2D_tmp = SNR_density_SG(:,:,i2);
    SNR_2D_tmp = SNR_SG(:,:,i2);
    
    [Ux_fit_SG(i2),Uy_fit_SG(i2),~,lin_ind] = retrieve_best_fit_from_SNR_2D(SNR_density_2D_tmp,Ux_SG_2D ,Uy_SG_2D ,0,0);
    SNR_SG_max(i2) = SNR_2D_tmp(lin_ind);
    SNR_density_SG_max(i2) = SNR_density_2D_tmp(lin_ind);
    
%     %get SNR dispersion in the 2D fit
%     SNR_SG_dispersion(i2) = get_dispersion_on_2D_fit(SNR_2D_tmp,Ux_SG_2D,Uy_SG_2D,lin_ind,dispersion_conf_area_mps);
%     SNR_density_SG_dispersion(i2) = get_dispersion_on_2D_fit(SNR_density_2D_tmp,Ux_SG_2D,Uy_SG_2D,lin_ind,dispersion_conf_area_mps);
    
    % display first guess result
    disp(['Ux_SG: ' num2str(Ux_fit_SG(i2)) ' m/s   Uy_SG: ' num2str(Uy_fit_SG(i2)) ' m/s   SNR density: ' ...
        num2str(SNR_density_SG_max(i2)) ' step K: ' num2str(K_steps_AV(i2))]);
    
end

% figure;
% subplot(1,3,1);
% plot(K_steps_AV,SNR_density_SG_max)
% subplot(1,3,2);
% plot(K_steps_AV,SNR_density_SG_dispersion)
% subplot(1,3,3);
% plot(SNR_density_SG_max,SNR_density_SG_dispersion,'.r')



% create second guess sctructure
SG_fit =  struct('signal_3D',signal_SG_3D,'noise_2D',noise_SG_3D,...
       'signal_nvalues_3D',signal_nvalues_3D,'noise_nvalues_3D',noise_nvalues_3D,...
       'SNR_3D',SNR_SG,'SNR_density_3D',SNR_density_SG,...
       'Ux_fit_Kstep',Ux_fit_SG,'Uy_fit_Kstep',Uy_fit_SG,...
       'SNR_max_Kstep',SNR_SG_max,'SNR_density_max_Kstep',SNR_density_SG_max,...
       'Ux_2D',Ux_SG_2D,'Uy_2D',Uy_SG_2D,'K_steps_AV',K_steps_AV,...
       'K_steps_minmax',fit_param.K_steps);
   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  generate output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
fit_out = struct('FG_fit',FG_fit,'SG_fit',SG_fit,'window_used',window_used);

end

