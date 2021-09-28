function [SNR,P_k,G] = nsp_doppler_shift_extraction(Spectrum,params,kval,U,V)

%This function calcalates the signal-to-noise ratio of the defined spectral
%signal function weighted by a function defining the linear dispersion of
%waves at constant wavenumber as a function of angle in  the presence of a current. 
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

%params: structure containing various relevant physical parameters, with
%the following fields (all units must be consistent):
%params.g - gravitational acceleration
%params.T - surface tension coefficient / density
%params.h - water depth
%params.omegaWidth - frequency width of weighting function (half width
%1/e^2)

%kval - wavenumber value

%U - current speed along the x-direction (1st dimension of WS.Pspec)

%V - current speed along the y-direction (2nd dimension of WS.Pspec)

%SpecCylinderCS (optional structure) - wave spectrum signal on a cylindrical surface of constant
%wavenumber. If provided as an input, this structure OVERIDES 'Spectrum'
%SpecCylinderCS.P_k: (omega, theta) wave spectrum signal on cylindrical surface of constant wavenumber
%SpecCylinderCS.thetaM: mesh of angle values [rad] on which spectrum is defined
%SpecCylinderCS.omegaM: mesh of frequency values [rad/s] on which spectrum is defined

%include2ndHarm (optional) - whether to include the 2nd harmonic in P_k. Set to FALSE by default.

%logFlag (optional) - whether to define P_k as the log of the spectrum. Set to TRUE by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OUTPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SNR - signal-to-noise ratio of the weighted spectral signal

%P_k - spectral signal on cylindrical surface of constant wavenumber

%thetaM - mesh of angle values on which P_k is defined (second dimension)

%omegaM - mesh of frequency values on which P_k is defined (first dimension)

%G - spectral weighting function.
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Define wave dispersion relation.
omegaFun = @(kx,ky) sqrt((params.g*sqrt(kx.^2 + ky.^2) + params.T*sqrt(kx.^2 + ky.^2).^3).*tanh(sqrt(kx.^2 + ky.^2)*params.h)) + U.*kx + V.*ky;


%Frequency width of weighting function (1/e^2 halfwidth)
a = params.omegaWidth;

if isfield(params,'kWidth')
    a_k = params.kWidth;
else
    a_k = max(Spectrum.dKx,Spectrum.dKy) * 2;
end

inds = ~isnan(Spectrum.power_Spectrum);%Discard indices where the spectrum contains NaN's

% P_k =    sqrt(Spectrum.power_Spectrum);%: (Kx,Ky,W) power Spectrum, i.e. x-direction is 1st dimension of power_Spectrum, y-direction is 2nd dimension, W (omega) is 3rd dimension.    
% KX =     Spectrum.Kx_3D;%: 3D Kx grid corresponding to Spectrum.power_Spectrum [rad/m]
% KY =     Spectrum.Ky_3D;%: 3D Ky grid corresponding to Spectrum.power_Spectrum [rad/m]
% W =      Spectrum.W_3D;%: 3D W grid corresponding to Spectrum.power_Spectrum [rad/sec]

P_k =    sqrt(Spectrum.power_Spectrum(inds));%: (Kx,Ky,W) power Spectrum, i.e. x-direction is 1st dimension of power_Spectrum, y-direction is 2nd dimension, W (omega) is 3rd dimension.    
KX =     Spectrum.Kx_3D(inds);%: 3D Kx grid corresponding to Spectrum.power_Spectrum [rad/m]
KY =     Spectrum.Ky_3D(inds);%: 3D Ky grid corresponding to Spectrum.power_Spectrum [rad/m]
W =      Spectrum.W_3D(inds);%: 3D W grid corresponding to Spectrum.power_Spectrum [rad/sec]


K = sqrt(KX.^2 + KY.^2);
order = 2;


if ~isnan(kval)
    P_k = P_k.*exp( -2*((K-kval)/a_k).^order);
end

%Define weighting function G
G1 = exp( -2*((W - omegaFun( KX, KY))/a).^order);
G2 = exp( -2*((W + omegaFun(-KX,-KY))/a).^order);
G = G1+G2;

P_k(~isfinite(P_k)) = 0;
InP = P_k.*G;

signal = sum(InP(:))/sum(G(:));
noise = sum(P_k(:).*(1-G(:)))/sum(1-G(:));
SNR = signal./noise;

end