function Scyl = cylinder_cross_section(Spectrum,dtheta,order)


if ~exist('order','var')
    order = 4;
end

Kx_2D = squeeze(Spectrum.Kx_3D(:,:,1));
Ky_2D = squeeze(Spectrum.Ky_3D(:,:,1));
W_1D = squeeze(Spectrum.W_3D(1,1,:));

theta_2D = atan2(-Ky_2D,-Kx_2D);
theta_3D = repmat(theta_2D,1,1,numel(W_1D));

thetaVals = -pi:dtheta:pi;

S_cyl = zeros(numel(W_1D),numel(thetaVals));

for ii = 1:numel(thetaVals)
    
    thetaDiff = angle(exp(1.i*theta_3D).*exp(-1.i*thetaVals(ii)));
    Theta_filt = exp(-2*(thetaDiff/dtheta).^order);
    Sfilt = Spectrum.power_Spectrum.*Theta_filt;    
    S_cyl(:,ii) = squeeze(nansum(Sfilt,[1,2]));

end

[thetaM,omegaM] = meshgrid(thetaVals,W_1D);
Scyl.P_k = S_cyl;
Scyl.thetaM = thetaM;
Scyl.omegaM = omegaM;


end