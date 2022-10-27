function out = find_current_depth_profile(wavenumbers,Ux,Uy,waterDepth)

%This function finds the depth profile of the currents, using a number of
%various methods (could be added to in the future). All methods work from 
%{k,Ux,Uy} triplets obtained from a 3D wave spectrum. So far, most methods
%only work in infinite depth. 
%To date, the methods used are:

% 1) Effective depth method (EDM): depths Z_eff are assigned to each wavenumber value
%based on an assumption to the functional form of the current profile.
%Three options are performed:

    %Linear profile: Z_eff = -1/(2*k)

    %Logarithmic profile: Z_eff = -1/(3.56*k)

    %Phillips wind spectrum Stokes drift profile: Z_eff = -pi/(8*k)
    %See e.g. Ø. Breivik et al. 2016 A Stokes drift approximation based on the
    %Phillips spectrum. Ocean Modelling 100, 49-56.

%For each case above, the mapped velocities are used to approximate an
%appropriate depth profile based on the assumed functional form.


%2) Polynomial effective depth method (PEDM): The method is an extension of
%the EDM method assuming a linear profile in the mapping. The method then
%corrects the velocities at each depth assuming a polynomial form to the
%current profile. 
%See: 
%Smeltzer, B. K., Æsøy, E., Ådnøy, A.,& Ellingsen, S. Å. (2019). An improved
%method for determining near-surface currents from wave dispersion measurements. 
%Journal of Geophysical Research: Oceans, 124. https://doi.org/10.1029/2019JC015202

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%wavenumbers - values at which Doppler shift currents are measured [rad/m]
%Ux - doppler shift velocity in the x-direction [m/s]
%Uy - doppler shift velocity in the y-direction [m/s]

%waterDepth (optional): set to Inf by default. Only used for the EDM-linear
%and PEDM methods above

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OUTPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The function rotates the horizontal coordinate system to 'along-flow'
%and 'cross-flow' directions instead of 'x' and 'y' the output current
%profile are provided in this coordinate system, with rotation angle 'phi'
%to transform back to the input coordinate system.

%out.global.wavenumbers - wavenumbers used in the method, same as input unless some of the Ux, Uy input values are NaN's
%out.global.U1 - Doppler shifts in along-flow direction
%out.global.U2 - Doppler shifts in cross-flow direction
%out.phi - rotation angle to transform back input coordinate system:
    %Ux = U1(inds)*cos(phi) - U2(inds)*sin(phi);
    %Uy = U1(inds)*sin(phi) + U2(inds)*cos(phi);


%out.EDM - structure with following substructures: 'lin','log','Phillips'
%All have the following fields:
%Z_eff: effective depths
%U1_fun: depth profile function handle in along-flow direction (depths are negative values)
%U2_fun: depth profile function handle in cross-flow direction (depths are negative values)
%eps1: RMS difference between forward calculated Doppler shifts (using U1_fun) and U1 values
%eps2: RMS difference between forward calculated Doppler shifts (using U2_fun) and U2 values

%out.EDM.lin.Usurf1 - surface velocity based on linear profile fit
%out.EDM.lin.Shear1 - constant vertical shear (vorticity) based on linear profile fit

%out.EDM.log.uStar1 - friction velocity based on log fit
%out.EDM.Phillips.uWind1 - wind velocity based on Phillips spectrum profile Stokes drift profile fit

%The above fields out.EDM.***.***2 reflect the same quantities in the cross flow direction 

%out.PEDM.Z_eff: effective depths
%out.PEDM.U1_fun: depth profile function handle in along-flow direction (depths are negative values)
%out.PEDM.U2_fun: depth profile function handle in cross-flow direction (depths are negative values)
%out.PEDM.eps1: RMS difference between forward calculated Doppler shifts (using U1_fun) and U1 values
%out.PEDM.eps2: RMS difference between forward calculated Doppler shifts (using U2_fun) and U2 values
%out.PEDM.verbose1 - infomation about the PEDM fits in along-flow direction - see pedm.m for more information
%out.PEDM.verbose2 - infomation about the PEDM fits in cross-flow direction - see pedm.m for more information

%If user has not provided information about the waterDepth, set to infinity
if ~exist('waterDepth','var')
    waterDepth = Inf;
end

if isempty(waterDepth)
    waterDepth = Inf;
end


indsX = ~isnan(Ux);
indsY = ~isnan(Uy);
inds = indsX & indsY;
wavenumbers = wavenumbers(inds);


U0 = mean(Ux(inds));
V0 = mean(Uy(inds));
%phi = atan2(V0,U0);
phi = 0;

U1 = Ux(inds)*cos(-phi) - Uy(inds)*sin(-phi);
U2 = Ux(inds)*sin(-phi) + Uy(inds)*cos(-phi);

out.global.wavenumbers = wavenumbers;
out.global.U1 = U1;
out.global.U2 = U2;
out.global.phi = phi;
out.global.Ux = Ux(inds);
out.global.Uy = Uy(inds);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%EDM%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z_eff_lin = -(2*wavenumbers).^-1.*tanh(abs(waterDepth)*wavenumbers);%Effective depths assuming linear profile.
Z_eff_log = -(3.56*wavenumbers).^-1;%Effective depths assuming logarithmic profile.
%Z_eff_Phillips = -pi*(8*wavenumbers).^-1;%Effective depths assuming Phillips Stokes drift profile.

out.EDM.lin.Z_eff = Z_eff_lin;
out.EDM.log.Z_eff = Z_eff_log;
%out.EDM.Phillips.Z_eff = Z_eff_Phillips;

% alpha = 0.008;
% g = 9.81;
%U_Phill = @(Uw,z) 2*alpha*Uw.*(exp(2*g./Uw.^2.*z) - ...
%    sqrt(-2*pi*g./Uw.^2.*z).*erfc(sqrt(-2*g./Uw.^2.*z)));%Phillips Stokes drift profile
kappa = 0.41;%Von Karman constant.

for i = 1:2
    switch i
        case 1
            uvals = U1;
        case 2
            uvals = U2;
    end

pLog = polyfit(log(wavenumbers),uvals,1);
pLin = polyfit(Z_eff_lin,uvals,1);
%Uwind = lsqcurvefit(U_Phill,10,Z_eff_Phillips,uvals);

usKap = pLog(1);
U0lnz0 = pLog(2) - usKap*log(2*1.78);

%Forward-calculated Doppler shifts.
uTildeLog = polyval(pLog,log(wavenumbers));
uTildeLin = -pLin(1)./(2*wavenumbers).*tanh(wavenumbers*waterDepth) + pLin(2);
%uTildePhill = 2*alpha*Uwind*(1 + sqrt(g./Uwind.^2./wavenumbers).*...
%    (atan(sqrt(g./Uwind.^2./wavenumbers)) - pi/2));

epsLog = rms(uTildeLog-uvals);
epsLin = rms(uTildeLin-uvals);
%epsPhill = rms(uTildePhill-uvals);

UlogFun = @(z) U0lnz0 - usKap*log(abs(z));
uStar = pLog(1)*kappa;

switch i
    case 1
        
        out.EDM.lin.Usurf1 = pLin(2);
        out.EDM.lin.Shear1 = pLin(1);
        out.EDM.lin.U1_fun = @(z) polyval(pLin,z);
        out.EDM.lin.eps1 = epsLin;
               
        out.EDM.log.U1_fun = UlogFun;
        out.EDM.log.eps1 = epsLog;
        out.EDM.log.uStar1 = uStar;

%         out.EDM.Phillips.U1_fun = @(z) U_Phill(Uwind,z);
%         out.EDM.Phillips.eps1 = epsPhill;
%         out.EDM.Phillips.uWind1 = Uwind;
        
    case 2
        out.EDM.lin.Usurf2 = pLin(2);
        out.EDM.lin.Shear2 = pLin(1);
        out.EDM.lin.U2_fun = @(z) polyval(pLin,z);
        out.EDM.lin.eps2 = epsLin;
       
        out.EDM.log.U2_fun = UlogFun;
        out.EDM.log.eps2 = epsLog;
        out.EDM.log.uStar2 = uStar;
% 
%         out.EDM.Phillips.U2_fun = @(z) U_Phill(Uwind,z);
%         out.EDM.Phillips.eps2 = epsPhill;
%         out.EDM.Phillips.uWind2 = Uwind;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PEDM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[pPEDM_1,pEDM_1,eps_PEDM_1,eps_EDM_1,verbose_1] = pedm(wavenumbers,U1,[],[],[],waterDepth);
[pPEDM_2,pEDM_2,eps_PEDM_2,eps_EDM_2,verbose_2] = pedm(wavenumbers,U2,[],[],[],waterDepth);

out.PEDM.U1_fun = @(z) polyval(pPEDM_1,z);
out.PEDM.U2_fun = @(z) polyval(pPEDM_2,z);
out.PEDM.eps1 = eps_PEDM_1;
out.PEDM.eps2 = eps_PEDM_2;
out.PEDM.verbose1 = verbose_1;
out.PEDM.verbose2 = verbose_2;

% out.EDM.lin.U1_fun = @(z) polyval(pEDM_1,z);
% out.EDM.lin.U2_fun = @(z) polyval(pEDM_2,z);
% out.EDM.lin.eps1 = eps_EDM_1;
% out.EDM.lin.eps2 = eps_EDM_2;



end