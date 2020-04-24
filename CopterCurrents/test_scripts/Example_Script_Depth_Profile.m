%Example script demonstrating inversion methods for finding the current
%depth profile from wave Doppler shift measurements.
%Infinite depth is assumed.

wavelength_min = 15.0;
wavelength_max = 100;
nWavenumbers = 20;%Number of wavenumbers

wavenumbers = linspace(2*pi/wavelength_max,2*pi/wavelength_min,nWavenumbers);

Uprof = 'exp';%'lin','log','exp','Phillips'...
noise = 1.0e-3;%Rad/s

switch Uprof
    case 'lin'
        U0 = 1.0;%Surface current
        S = 0.2;%Shear strength
        util = -S./(2*wavenumbers) + U0;
        Ufun = @(z) U0 + S*z;
        
    case 'log'
        U0 = 0.4;%Surface velocity.
        z0 = 1e-2;%Roughness length
        uStar = 0.02;%Friction velocity
        kappa = 0.4;
        util = U0 - uStar/kappa*log(1./(2*wavenumbers*1.78)./z0);
        Ufun = @(z) U0 - uStar/kappa*log(abs(z)./z0);
       
    case 'exp'
        U0 = 0.25;
        z1 = 2.0;
        Ufun = @(z) U0*exp(z/z1);
        util = 2*U0*wavenumbers./(1/z1 + 2*wavenumbers);
    case 'Phillips'
        alpha = 0.008;%Constant
        g = 9.81;
        Uw = 20;%Wind speed 
        Ufun = @(z) 2*alpha*Uw.*(exp(2*g./Uw.^2.*z) - ...
            sqrt(-2*pi*g./Uw.^2.*z).*erfc(sqrt(-2*g./Uw.^2.*z)));
        
        util = 2*alpha*Uw*(1 + sqrt(g./Uw.^2./wavenumbers).*...
    (atan(sqrt(g./Uw.^2./wavenumbers)) - pi/2));
        
        
end

%Add the noise.
util = util + randn(1,nWavenumbers)*noise./wavenumbers;


%uTildeLog = polyval(pLog,log(wavenumbers));
%uTildeLin = -pLin(1)./(2*wavenumbers).*tanh(wavenumbers*waterDepth) + pLin(2);
%uTildePhill = 2*alpha*Uwind*(1 + sqrt(g./Uwind.^2./wavenumbers).*...
%    (atan(sqrt(g./Uwind.^2./wavenumbers)) - pi/2));

Z_eff_lin = -(2*wavenumbers).^-1;%Effective depths assuming linear profile.
Z_eff_log = -(3.56*wavenumbers).^-1;%Effective depths assuming logarithmic profile.
Z_eff_Phillips = -pi*(8*wavenumbers).^-1;%Effective depths assuming Phillips Stokes drift profile.

zvect = linspace(Z_eff_lin(1),0,400);

figure(1);subplot(1,2,1);plot(util,wavenumbers,'sk');xlabel('[m/s]');ylabel('k [rad/m]');
figure(1);subplot(1,2,2);plot(...
    util,Z_eff_lin,'o',...
    util,Z_eff_log,'d',...
    Ufun(zvect),zvect,'--k');
xlabel('[m/s]');ylabel('Depth [m]');
legend('Linear','Log','True Profile','Location','se');
drawnow;

%Now find the depth profile.
out = find_current_depth_profile(wavenumbers,util,0*util);
plot_depth_profile_results(out,3,Ufun);







