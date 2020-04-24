function [SNR,SpecCylinderCS,G] = k_cylinder_nsp(Spectrum,params,kval,U,V,SpecCylinderCS,include2ndHarm,logFlag)

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

%Set defaults if inputs are not provided.
if isempty(include2ndHarm)
    include2ndHarm = 0;
end

if isempty(SpecCylinderCS)
   P_k = []; 
else
    P_k = SpecCylinderCS.P_k;
    thetaM = SpecCylinderCS.thetaM;
    omegaM = SpecCylinderCS.omegaM;
end

if isempty(logFlag)
   logFlag = 0; 
end


%If the user has not provided 'SpecCylinderCS' as input, calculate what's required from 'Spectrum'.
if isempty(P_k)

[~,iwm] = max(Spectrum.W_3D(1,1,:));
[thetaM,omegaiM] = meshgrid(linspace(0,2*pi,400),(1:iwm));

omegaM = (omegaiM-1)*Spectrum.dW;%Omega-vals in rad/s (omegaiM are indices).

P_k = zeros(size(thetaM));
P_k2 = zeros(size(thetaM));%Second harmonic

%Spectrum.Kx_3D,Ky_3D,W_3D may have defined by 'fftshift()', or defined in
%the MATLAB convention of 'fftn()'. We handle both cases for robustness,
%and define a 'shifting' function.
if abs(Spectrum.Kx_3D(1,1,1))>0
    fftshiftFun = @(A) A;%If Spectrum has already been fftshifted, do nothing
else
    fftshiftFun = @(A) fftshift(A);%Otherwise, fftshift() the fields.
end

KX = squeeze(Spectrum.Kx_3D(:,:,1));%2D mesh
KY = squeeze(Spectrum.Ky_3D(:,:,1));%2D mesh
%[KX,KY] = meshgrid(WS.kxv,WS.kyv);

for i = 1:size(omegaiM,1)
    
    if logFlag
        
        P_ki = 10*log10(squeeze(Spectrum.power_Spectrum(:,:,omegaiM(i))));
        ind2 = (omegaiM(i)-1)*2 +1;%if ind2 > size(WS.Pspec,3);ind2 = 1;end
        %Don't include 2nd harmonic in the case of aliasing (Handling for
        %aliasing could be included in the future)
        if ind2 < size(Spectrum.power_Spectrum,3);P_ki2 = 10*log10(squeeze(Spectrum.power_Spectrum(:,:,ind2)));end
    else
        %
        P_ki = (squeeze(Spectrum.power_Spectrum(:,:,omegaiM(i))));
        ind2 = (omegaiM(i)-1)*2 +1;%if ind2 > size(WS.Pspec,3);ind2 = 1;end
        %P_ki2 = (squeeze(WS.Pspec(:,:,ind2)));
        %Same as above w.r.t. aliasing.
        if ind2 < size(Spectrum.power_Spectrum,3);P_ki2 = 10*(squeeze(Spectrum.power_Spectrum(:,:,ind2)));end

    end
    
    P_k(i,:) = interp2(fftshiftFun(KY),fftshiftFun(KX),fftshiftFun(P_ki),kval*sin(thetaM(1,:)),kval*cos(thetaM(1,:)));
    if include2ndHarm
        P_k2(i,:) = interp2(fftshiftFun(KY),fftshiftFun(KX),fftshiftFun(P_ki2),2*kval*cos(thetaM(1,:)),2*kval*sin(thetaM(1,:)));
    end
end


if include2ndHarm
P_k2(P_k2>P_k) = min(P_k2(isfinite(P_k2(:))));%Nix components of the 2nd harmonic which are greater than the first harmonic component.
try P_k2(isnan(P_k2)) = min(P_k2(isfinite(P_k2(:))));catch;P_k2 = 0;end %Also nix NaN components. If everything is NaN, set P_k2 = 0;
P_k = P_k + P_k2;
end

if logFlag
    P_k = P_k - min(P_k(isfinite(P_k)));
end


end

%Define wave dispersion relation.
omegaFun_kth = @(k,theta) sqrt((params.g*k + params.T*k.^3).*tanh(k*params.h)) + U*cos(theta).*k + V*sin(theta).*k;


%Frequency width of weighting function (1/e^2 halfwidth)
a = params.omegaWidth;

%Define weighting function G
G1 = exp( -2*((omegaM - omegaFun_kth( kval,thetaM))/a).^2);
G2 = exp( -2*((omegaM + omegaFun_kth(-kval,thetaM))/a).^2);
G = G1+G2;

P_k(~isfinite(P_k)) = 0;
InP = P_k.*G;


signal = sum(InP(:))/sum(G(:));
noise = sum(P_k(:).*(1-G(:)))/sum(1-G(:));
if logFlag
    SNR = 10^((signal-noise)/10);
else
    SNR = signal./noise;
end

SpecCylinderCS.P_k = P_k;
SpecCylinderCS.thetaM = thetaM;
SpecCylinderCS.omegaM = omegaM;


end