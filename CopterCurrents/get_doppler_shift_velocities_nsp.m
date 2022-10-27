function out_fit = get_doppler_shift_velocities_nsp(Spectrum,fit_param,Properties,verboseFig)


%This function extracts Doppler shift velocities representing the
%alteration in wave phase velocity due to a depth-dependent current.
%The function is based on the normalized scalar product method, described in the following
%article:

%Smeltzer, B. K., Æsøy, E., Ådnøy, A.,& Ellingsen, S. Å. (2019). An improved
%method for determining near-surface currents from wave dispersion measurements. 
%Journal of Geophysical Research: Oceans, 124. https://doi.org/10.1029/2019JC015202

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Spectrum structure with the the 3D wave spectrum (across horizontal space
% and time dimensions) containing the following fields:
%    Spectrum.power_Spectrum: (Kx,Ky,W) power Spectrum, i.e. x-direction is 1st dimension of power_Spectrum, y-direction is 2nd dimension, W (omega) is 3rd dimension.    
%    Spectrum.Kx_3D: 3D Kx grid corresponding to Spectrum.power_Spectrum [rad/m]
%    Spectrum.Ky_3D: 3D Ky grid corresponding to Spectrum.power_Spectrum [rad/m]
%    Spectrum.W_3D: 3D W grid corresponding to Spectrum.power_Spectrum [rad/sec]
%    Spectrum.dKx: Kx resolution [rad/m]
%    Spectrum.dKy: Ky resolution [rad/m]
%    Spectrum.dW: W resolution [rad/sec]

%fit_param.Ux_2D: 2D mesh with U_x values
%fit_param.Uy_2D: 2D mesh with U_y values
%fit_param.wavenumbers: list of wavenumbers in [rad/length] at which Doppler shift velocities will be extracted
%fit_param.include2ndHarmonic (optional): whether to include 2nd harmonic of the spectrum in the fit (set to false by default)
%fit_param.logFlag (optional): whether to do the fit in log space (false by default)
%fit_param.omegaWidthFun: function handle as a function of wavenumber i.e.
%@(k) f(k)...., specifying frequency width of the weighting function in
%frequency-angle space (constant wavenumber). Width is half-width 1/e^2
%point of a Gaussian function

%The following parameters involve post-processing of the Doppler shifts:
%fit_param.SNR_filter: whether to use a signal-to-noise filter (false by
%fit_param.SNR_threshold: threshold signal-to-noise value for above filter (set to 2.0 by default)
%fit_param.Peak_filter: whether to use a multiple peaks filter (false by default)
%fit_param.Peak_threshold: peak threshold of maximum value (0.5 by default)
%fit_param.Outlier_filter: whether to use an outlier filter (quartile-based) (false by default)

%Properties.g - gravitational acceleration
%Properties.T - surface tension coefficient / density
%Properties.h - water depth
%NOTE: All units must be consistent in dimensions.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OUTPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%out_fit.wavenumbers: list of wavenumbers (same as fit_param.wavenumbers)
%out_fit.Ux: list of raw Ux fit values
%out_fit.Uy: list of raw Uy fit values
%out_fit.SNR_max: list peak SNR values at each wavenumber
%out_fit.Ux_2D: same as fit_params.Ux_2D
%out_fit.Uy_2D: same as fit_params.Uy_2D
%out_fit.Ux_filt: filtered Ux fit values
%out_fit.Uy_filt: filtered Uy fit values

%The following is an array 'out_fit.verbose' containing further information for each wavenumber
%concerning the filtering.
%out_fit.verbose.k : wavenumber
%out_fit.verbose.SNR_2D: SNR values on mesh defined by out_fit.Ux_2D, out_fit_Uy_2D
%out_fit.verbose.Peaks_2D: logical array identifying peaks

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parse inputs, set defaults if needed.

if ~exist('verboseFig','var')
    verboseFig = 0;
end


if ~isfield(fit_param,'include2ndHarmonic')
   fit_param.include2ndHarmonic = 0; 
end

if ~isfield(fit_param,'logFlag')
   fit_param.logFlag = 0; 
end

if ~isfield(fit_param,'SNR_filter')
   fit_param.SNR_filter = 0; 
end

if ~isfield(fit_param,'SNR_threshold')
   fit_param.SNR_threshold = 2.0; 
end

if ~isfield(fit_param,'Peak_filter')
   fit_param.Peak_filter = 0; 
end

if ~isfield(fit_param,'Peak_threshold')
   fit_param.Peak_threshold = 0.5; 
end

if ~isfield(fit_param,'Outlier_filter')
   fit_param.Outlier_filter = 0; 
end


%Initial the output structure:
out_fit.wavenumbers = fit_param.wavenumbers;
out_fit.Ux = zeros(1,numel(fit_param.wavenumbers));
out_fit.Uy = zeros(1,numel(fit_param.wavenumbers));
out_fit.SNR_max = zeros(1,numel(fit_param.wavenumbers));
out_fit.Ux_2D = fit_param.Ux_2D;
out_fit.Uy_2D = fit_param.Uy_2D;
out_fit.verbose = [];

h = Properties.h;
g = Properties.g;
T = Properties.T;
omegaWidthFun = fit_param.omegaWidthFun;
wavenumbers = fit_param.wavenumbers;

%First, come up with an initial guess, based on the current over all
%wavenumbers.

% P = struct('h',h,'g',g,'T',T,'omegaWidth',omegaWidthFun(mean(wavenumbers)),'kWidth',0.05);
% 
% %Evaluate the NSP on grid of points to avoid other local maxima
% snrG = zeros(size(fit_param.Ux_2D));
% for i = 1:size(fit_param.Ux_2D,1)
%     for j = 1:size(fit_param.Ux_2D,2)
%         [snr_ij,~,~] = nsp_doppler_shift_extraction(Spectrum,P,NaN,fit_param.Ux_2D(i,j),fit_param.Uy_2D(i,j));
%         snrG(i,j) = snr_ij;
%     end
% end
% [~,im] = max(snrG(:));
% 
% figure(1002);
% %subplot(1,3,1);
% imagesc(fit_param.Uy_2D(1,:),fit_param.Ux_2D(:,1),snrG);colorbar;axis image;
% xlabel('U_y [m/s]');ylabel('U_x [m/s]');%title(sprintf('SNR: k = %.3f rad/m',wavenumbers(jj)));
% drawnow;
% 
% cDoppGuess = [fit_param.Ux_2D(im),fit_param.Uy_2D(im)]

for jj = 1:numel(fit_param.wavenumbers)
  
P = struct('h',h,'g',g,'T',T,'omegaWidth',omegaWidthFun(wavenumbers(jj)),'kWidth',fit_param.kWidth,'logFlag',fit_param.logFlag);

snrG = zeros(size(fit_param.Ux_2D));
for i = 1:size(fit_param.Ux_2D,1)
    for j = 1:size(fit_param.Ux_2D,2)
        [snr_ij,~,~] = nsp_doppler_shift_extraction(Spectrum,P,wavenumbers(jj),fit_param.Ux_2D(i,j),fit_param.Uy_2D(i,j));
        snrG(i,j) = snr_ij;
    end
end
[~,im] = max(snrG(:));

if verboseFig
figure(1002);
%subplot(1,3,1);
%imagesc(fit_param.Uy_2D(1,:),fit_param.Ux_2D(:,1),snrG);colorbar;axis image;
contourf(-fit_param.Uy_2D,-fit_param.Ux_2D,snrG,20);colorbar;axis image;
xlabel('$U_y$ [m/s]','Interpreter','latex');ylabel('$U_x$ [m/s]','Interpreter','latex');
title(sprintf('SNR: $k = %.2f$ rad/m',wavenumbers(jj)),'Interpreter','latex','FontWeight','normal');
%title(sprintf('SNR: wavenumber = %.2f rad/m',wavenumbers(jj)),'FontWeight','normal');
grid on;
drawnow;
end

cDoppGuess = [fit_param.Ux_2D(im),fit_param.Uy_2D(im)];

options = optimset('Display','off');%Turn off display output
Ufit = fminsearch(@(U) 1-nsp_doppler_shift_extraction(Spectrum,P,wavenumbers(jj),U(1),U(2)),cDoppGuess,options);

[SNR_fit,P_k,G] = nsp_doppler_shift_extraction(Spectrum,P,wavenumbers(jj),Ufit(1),Ufit(2));

Ux(jj) = Ufit(1);
Uy(jj) = Ufit(2);
SNR_max(jj) = SNR_fit;


verbose(jj) = struct('k',fit_param.wavenumbers(jj),'SNR_2D',snrG,'Peaks_2D',find_peaks(snrG));

end


%Update the output structure.
out_fit.Ux = Ux;
out_fit.Uy = Uy;
out_fit.SNR_max = SNR_max;
out_fit.verbose = verbose;


%%%%%% POST-PROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Filter the ouputs as requested


out_fit.Ux_filt = out_fit.Ux;
out_fit.Uy_filt = out_fit.Uy;
if fit_param.SNR_filter
    SNR_inds = out_fit.SNR_max > fit_param.SNR_threshold;
    out_fit.Ux_filt(~SNR_inds) = NaN;
    out_fit.Uy_filt(~SNR_inds) = NaN;
end



if fit_param.Peak_filter
    indsMultPeaks = peak_filter(out_fit.verbose,fit_param.Peak_threshold);
    out_fit.Ux_filt(indsMultPeaks) = NaN;
    out_fit.Uy_filt(indsMultPeaks) = NaN;   
end


if fit_param.Outlier_filter
    indsX_keep = quartile_filter(out_fit.wavenumbers,out_fit.Ux_filt);
    indsY_keep = quartile_filter(out_fit.wavenumbers,out_fit.Uy_filt);
    out_fit.Ux_filt(~indsX_keep) = NaN;
    out_fit.Uy_filt(~indsY_keep) = NaN;
end



end


%Peak filter. returns indices of structure array where there are multiple
%peaks in the SNR_2D field exceeding a threshold value of the max SNR-value
function indsMultiplePeaks = peak_filter(S,peakThresh)

indsMultiplePeaks = zeros(1,numel(S));

for i = 1:numel(S)
peakVals = S(i).SNR_2D(S(i).Peaks_2D);
indsMultiplePeaks(i) = numel(find(peakVals/max(peakVals) > peakThresh))>1;
end

indsMultiplePeaks = logical(indsMultiplePeaks);

end

%Filter based on quartiles
function indsOutlier = quartile_filter(x,y)

indsNaN = isnan(y);

x = x(~indsNaN);
y = y(~indsNaN);

%yFilt = y;
pct = polyfit(x,y,1);

delCtil = y-polyval(pct,x);
delCtil_srt = sort(delCtil);

Q1 = median(delCtil_srt(1:round(numel(delCtil_srt)/2)));Q3 = median(delCtil_srt(round(numel(delCtil_srt)/2):end));
IQR = Q3-Q1;
lowerFence = Q1-1.5*IQR;upperFence = Q3+1.5*IQR;

indsKeep = and(delCtil>lowerFence,delCtil<upperFence);

%yFilt(~indsKeep) = NaN;

indsOutlier = ones(size(indsNaN));
indsOutlier(~indsNaN) = indsKeep;

end




