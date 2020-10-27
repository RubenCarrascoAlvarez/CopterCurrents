function Error_mps_vec = Error_depth_fit_model_order_n(beta,x)
% radarradar_backscatter_model_3D defines a radar backscatter model 
%   VV or HH antena, including the wind.
%   
%   Ux = x(:,1);  
%   Uy = x(:,2);  
%   Depth = x(:,3);  
%
% model definiton
%
%     order 1
%         Ux_fit = beta(1) + (beta(2) .* depth);
%         Uy_fit = beta(3) + (beta(4) .* depth); 
%         
%     order 2
%         Ux_fit = beta(1) + (beta(2) .* depth) + (beta(3) .* depth.^2);
%         Uy_fit = beta(4) + (beta(5) .* depth) + (beta(6).*  depth.^2); 
%         
%     order 3    
%         Ux_fit = beta(1) + (beta(2) .* depth) + (beta(3) .* depth.^2) + (beta(4) .* depth.^3);
%         Uy_fit = beta(5) + (beta(6) .* depth) + (beta(7).*  depth.^2) + (beta(8) .* depth.^3); 
%         
%     order 4    
%         Ux_fit = beta(1) + (beta(2) .* depth) + (beta(3) .* depth.^2) + (beta(4) .* depth.^3) + (beta(5) .* depth.^4);
%         Uy_fit = beta(6) + (beta(7) .* depth) + (beta(8).*  depth.^2) + (beta(9) .* depth.^3) + (beta(10) .* depth.^4);  
%         
%     order 5    
%         Ux_fit = beta(1) + (beta(2) .* depth) + (beta(3) .* depth.^2) + (beta(4)  .* depth.^3) + (beta(5)  .* depth.^4) + (beta(6)  .* depth.^5);
%         Uy_fit = beta(7) + (beta(8) .* depth) + (beta(9).*  depth.^2) + (beta(10) .* depth.^3) + (beta(11) .* depth.^4) + (beta(12) .* depth.^5) ; 


Ux_measured = x(:,1);  
Uy_measured = x(:,2);  
depth = x(:,3);  

% compute beta coefifcients
[Ux_fit,Uy_fit] = compute_beta_coef_depth_fit(beta,depth);

% compute error
Error_Ux = Ux_fit - Ux_measured;
Error_Uy = Uy_fit - Uy_measured;
Error_mps_vec = sqrt(Error_Ux.^2 + Error_Uy.^2);

              
end



