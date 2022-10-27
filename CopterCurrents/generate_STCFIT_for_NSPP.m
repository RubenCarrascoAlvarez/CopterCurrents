function STCFIT = generate_STCFIT_for_NSPP(STCFIT,wavenumbers,include2ndHarmonic,logFlag,omegaWidthFun,...
    SNR_filter,SNR_threshold,Peak_filter,Peak_threshold,Outlier_filter)

%STCFIT: structure containing relevant parameters for the current fit.

%wavenumbers
%List of wavenumber values at which to extract Doppler shift velocities
%[rad/m]. If empty, values will be chosen between the values of
%STCFIT.fit_param.K_limits steps of the k pixel spacing of each window (assumed to be
%square)
if isempty(wavenumbers)
%wavenumbers = linspace(STCFIT.fit_param.K_limits(1),STCFIT.fit_param.K_limits(2),20);
wavenumbers = STCFIT.fit_param.K_limits(1):2*pi/STCFIT.Windows.sq_size_m:STCFIT.fit_param.K_limits(2);
end


% include2ndHarmonic
% (optional): whether to include 2nd harmonic of the spectrum in the fit (false by default)
if isempty(include2ndHarmonic)
    include2ndHarmonic = 0;
end

%logFlag
% (optional): whether to do the fit in log space (false by default)
if isempty(logFlag)
logFlag = 0;
end

% (optional) omegaWidthFun: function handle as a function of wavenumber i.e.
%@(k) f(k)...., specifying frequency width of the weighting function in
%frequency-angle space (constant wavenumber). Width is half-width 1/e^2
%point of a Gaussian function.
if isempty(omegaWidthFun)
    dOmega = 2*pi/((STCFIT.Generic.time_stamp(end)-STCFIT.Generic.time_stamp(1))*3600*24);
omegaWidthFun = @(k) dOmega + 0.1*k;
end

%The following OPTIONAL parameters involve post-processing of the Doppler shifts:
%SNR_filter: whether to use a signal-to-noise filter (false by default)
if isempty(SNR_filter)
    SNR_filter = 0;
end

%SNR_threshold: threshold signal-to-noise value for above filter (set to 2.0 by default)
if isempty(SNR_threshold)
    SNR_threshold = 2.0;
end

%Peak_filter: whether to use a multiple peaks filter (false by default)
if isempty(Peak_filter)
    Peak_filter = 0;
end

%Peak_threshold: peak threshold of maximum value (0.5 by default)
if isempty(Peak_threshold)
    Peak_threshold = 0.5;
end

%Outlier_filter: whether to use an outlier filter (quartile-based) (false by default)
if isempty(Outlier_filter)
    Outlier_filter = 0.0;
end


STCFIT.NSPP_fit_param.wavenumbers = wavenumbers;
STCFIT.NSPP_fit_param.include2ndHarmonic = include2ndHarmonic;
STCFIT.NSPP_fit_param.logFlag = logFlag;
STCFIT.NSPP_fit_param.omegaWidthFun = omegaWidthFun;
STCFIT.NSPP_fit_param.SNR_filter = SNR_filter;
STCFIT.NSPP_fit_param.SNR_threshold = SNR_threshold;
STCFIT.NSPP_fit_param.Peak_filter = Peak_filter;
STCFIT.NSPP_fit_param.Peak_threshold = Peak_threshold;
STCFIT.NSPP_fit_param.Outlier_filter = Outlier_filter;


end