function [pPEDM_out,pEDM_out,eps_PEDM_out,eps_EDM_out,verbose] = pedm(k_vect,c_til,nMax_vals,deltaz_T_vals,deltaz_B_vals,waterDepth)

%This function is an implementation of the Polynomial effective
%depth method (PEDM) developed in the article:

%Smeltzer, B. K., Æsøy, E., Ådnøy, A.,& Ellingsen, S. Å. (2019). An improved
%method for determining near-surface currents from wave dispersion measurements. 
%Journal of Geophysical Research: Oceans, 124. https://doi.org/10.1029/2019JC015202

%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Units of inputs must be consistent, given as quantity in [...]

%k_vect - list of wavenumber values [rad/length]
%c_til - list of Doppler shift velocities at wavenumners k_vect [length/time]

%PEDM parameters (all positive values):
%nMax_vals %Polynomial order [unity]
%deltaz_T_vals %Depth interval (shallow end) for extrapolation [length]
%deltaz_B_vals %Depth interval (deep end) for extrapolation [length]

%If multiple values are given for the above PEDM parameters, the function will
%loop over all combinations. The outputs returned will reflect the
%parameter combination giving a minimum value of 'eps_PEDM', the RMS
%difference between 'c_til' and the forward calculated Doppler shifts based
%on the PEDM current profile. If no values are provided, the following
%defaults will be used:
%nMax_vals: (2:nm) where 'nm' is  min(12,round(numel(k_vect)/2)) 
%deltaz_T_vals: linspace(0.01,0.2,20)*depthRange where 'depthRange' is the range of mapped depths  
%deltaz_B_vals: linspace(0.02,0.8,20)*depthRange


%If nMax, deltaz_T, and/or deltaz_B contain multiple values, all parameter
%combinations will be run, with the outputs reflecting the combination with
%the smallest value of eps_PEDM.

%waterDepth (optional): water depth [length]. If not provided or [], infinite
%depth is assumed. %NOTE: this implementation of the PEDM has not been
%exclusively tested for cases of finite water depth. Use with some caution
%if min(k_vect*waterDepth)<pi (though it is expected to perform reasonably
%accurately for moderate depths - see above article for further details).

%%%%%%%%%%%%%% OUTPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%pPEDM: polynomial coefficients of the PEDM profile (with minimum eps_PEDM)
%pEDM: polynomial coefficients of the EDM profile (with minimum eps_PEDM)
%eps_PEDM: [length/time] RMS difference between c_til and forward calculated shifts based on PEDM profile
%eps_EDM: [length/time] RMS difference between c_til and forward calculated shifts based on EDM profile
%verbose: structure containing information of all polynomial fits over all PEDM parameter combinations

% pedm is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% pedm is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% <http://www.gnu.org/licenses/>.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define Root-mean-square function for later use
rms = @(x) sqrt(mean(x.^2));

%If user has not provided information about the waterDepth, set to infinity
if ~exist('waterDepth','var')
    waterDepth = Inf;
end

if isempty(waterDepth)
    waterDepth = Inf;
end

%Turn off warning about polyfits (will turn on after completion of the
%function)
warning('off','MATLAB:polyfit:RepeatedPointsOrRescale');

%Handling for inadequate input data.
if isempty(k_vect) || isempty(c_til)
    pPEDM_out = NaN;
    pEDM_out = NaN;
    eps_PEDM_out = NaN;
    eps_EDM_out = NaN;
    verbose = NaN;
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now we start with the PEDM. The steps follow those in the manuscript in
%section 2.1.2 Effect of limitations of measured Doppler shifts
%Define the parameters.
%Calculate effective depths of Doppler shift velocities based on assumption
%of a linear profile

Z_eff = -(2*k_vect).^-1.*tanh(abs(waterDepth)*k_vect);


%Set default values for PEDM parameter combinations if inputs left blank
if isempty(nMax_vals)
    nm = min(12,round(numel(k_vect)/2));
    nMax_vals = (0:nm);
end
if isempty(deltaz_T_vals)
    depthRange = abs(Z_eff(1)-Z_eff(end));
    deltaz_T_vals = linspace(0.01,0.2,20)*depthRange;    
end
if isempty(deltaz_B_vals)
    depthRange = abs(Z_eff(1)-Z_eff(end));
    deltaz_B_vals = linspace(0.02,0.8,20)*depthRange;
end


z_c = max(4*min(Z_eff),-abs(waterDepth));%Cutoff depth chosen as 4 times the deepest mapped depth. (Set to water depth if depth if shallower) 
%z_c is NEGATIVE by convention here


%We loop over all PEDM parameter combinations. First initialize
%eps_PEDM_out as well as the other outputs.

eps_PEDM_out = Inf;
eps_EDM_out = Inf;

pPEDM_out = NaN;
pEDM_out = NaN;


verbose = [];
combo = 0;


for nMax = nMax_vals
    for deltaz_B = deltaz_B_vals
        for deltaz_T = deltaz_T_vals


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1: Fit the mapped Doppler shifts to a polynomial of order nMax.
p1 = polyfit(Z_eff,c_til,nMax);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 2: Create additional velocity-depth pairs by linearly extrapolating
%up to the surface and down to cutoff depth z_c

zB = linspace(0,deltaz_B,100);
zT = linspace(-deltaz_T,0,100);

pTop = polyfit(Z_eff(end)+zT,polyval(p1,Z_eff(end)+zT),1);%Linear fit within depth interval deltaz_T from top of range of mapped depths
pBottom = polyfit(Z_eff(1) + zB,polyval(p1,Z_eff(1) + zB),1);%Linear fit within depth interval deltaz_B from bottom of range of mapped depths

depthsExBtm = z_c:deltaz_B:Z_eff(1)-deltaz_B;%Depths below deepest mapped depths
depthsExTop = [Z_eff(end)+deltaz_T:deltaz_T:0,0];%Depth values above shallowest mapped depth

zEx = [depthsExBtm,Z_eff,depthsExTop];%Expanded range of depths, including the depths where the current is extrapolated
cTilEx = [polyval(pBottom,depthsExBtm),c_til,polyval(pTop,depthsExTop)];%Expanded range of current values, including extrapolated 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 3: Perform a second polynomial fit on the expanded set of points to
%produce the profile considered U_EDM.
pEDM = polyfit(zEx,cTilEx,nMax);

%Determine if profile is monotonic with range of mapped depths.
nv = nMax:-1:1;
pczEx_p = nv.*pEDM(1:end-1);%Derivative polynomial coefficients
rts = roots(pczEx_p);%Find roots
rts = rts(imag(rts)==0);%Discard complex roots

%If there roots of the derivative within ranges of mapped depths, don't go
%further.
if any(rts>Z_eff(1) & rts<0)
    continue
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 4: Scale polynomial coefficients defining U_EDM by n! as in equation (8) in the article.
pPEDM_i = pEDM./factorial(nMax:-1:0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 5: Create a new set of linearly extrapolated points down to z_c based 
%on the average shear of of the above polynomial function in a depth
%interval deltaz_B/2 at the deep end of the regime.

zB2 = linspace(0,deltaz_B/2,100);
pBottom2 = polyfit(Z_eff(1) + zB,polyval(pPEDM_i,Z_eff(1) + zB),1);%Linear fit within depth interval deltaz_B from bottom of range of mapped depths

depthsExBtm2 = z_c:deltaz_B/2:Z_eff(1)-deltaz_B;%Depths below deepest mapped depths

zEx2 = [depthsExBtm2,Z_eff,depthsExTop];%Expanded range of depths, including the depths where the current is extrapolated
Uvals = [polyval(pBottom2,depthsExBtm2),polyval(pPEDM_i,Z_eff),polyval(pPEDM_i,depthsExTop)];%Expanded range of current values, including extrapolated 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 6: Perform a final polynomial fit on the expanded set of points.
pPEDM = polyfit(zEx2,Uvals,nMax);%FINAL POLYNOMIAL COEFFICIENTS (OUTPUT)

% Calculate Doppler shifts assuming U_EDM or U_PEDM as the current profile, using the forward problem. 
%See equation (9) and text above for explanation.

%Calculate c_til using equation (5) summation, or using imcomplete gamma
%functions, to only integrate down to depth z_c.
c_tilEDM = 0;
c_tilPEDM = 0;

if nMax<2

for n = 0:nMax
%Use equation (5).
c_tilEDM = c_tilEDM + factorial(n).*pEDM(end-n).*(-1./(2*k_vect)).^n;
c_tilPEDM = c_tilPEDM + factorial(n).*pPEDM(end-n).*(-1./(2*k_vect)).^n;

end

else

for n = 0:nMax
%Use incomplete gamma function
 %c_tilEDM = c_tilEDM + factorial(n).*pEDM(end-n).*(-1./(2*k_vect)).^n;
 c_tilEDM = c_tilEDM + (-1/2)^n*k_vect.^(-n)*pEDM(end-n).*(gammainc(-2*k_vect*z_c,1+n,'lower'))*gamma(n+1);
 %c_tilPEDM = c_tilPEDM + factorial(n).*pPEDM(end-n).*(-1./(2*k_vect)).^n;
 c_tilPEDM = c_tilPEDM + (-1/2)^n*k_vect.^(-n)*pPEDM(end-n).*(gammainc(-2*k_vect*z_c,1+n,'lower'))*gamma(n+1);

end
end

%Calculate RMS differences (equation (9) in mansuscript)
eps_EDM = rms(c_til-c_tilEDM);%RMS difference for U_EDM
eps_PEDM = rms(c_til-c_tilPEDM);%RMS difference for U_PEDM. 
%Paremeters nMax, deltaz_T, and deltaz_B are chosen to minimize eps_PEDM in
%practice.
if eps_PEDM<eps_PEDM_out
    eps_PEDM_out = eps_PEDM;
    %eps_EDM_out = eps_EDM;
    pPEDM_out = pPEDM;
    %pEDM_out = pEDM;
end

if eps_EDM<eps_EDM_out
    eps_EDM_out = eps_EDM;
    pEDM_out = pEDM;
end


combo = combo + 1;
verbose(combo).pPEDM = pPEDM;
verbose(combo).pEDM = pEDM;
verbose(combo).eps_PEDM = eps_PEDM;
verbose(combo).eps_EDM = eps_EDM;
verbose(combo).nMax = nMax;
verbose(combo).deltaz_B = deltaz_B;
verbose(combo).deltaz_T = deltaz_T;


        end
    end
end

warning('on','MATLAB:polyfit:RepeatedPointsOrRescale');

end