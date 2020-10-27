function [Ux_fit,Uy_fit] = compute_beta_coef_depth_fit(beta,depth)
% model definiton
% Ux_fit = bx1 + bx2*depth + bx3*depth^2 
% Uy_fit = by1 + by2*depth + by3*depth^2 


switch length(beta)
    
    case 4
        Ux_fit = beta(1) + (beta(2) .* depth);
        Uy_fit = beta(3) + (beta(4) .* depth); 
        
    case 6
        Ux_fit = beta(1) + (beta(2) .* depth) + (beta(3) .* depth.^2);
        Uy_fit = beta(4) + (beta(5) .* depth) + (beta(6).*  depth.^2); 
        
    case 8    
        Ux_fit = beta(1) + (beta(2) .* depth) + (beta(3) .* depth.^2) + (beta(4) .* depth.^3);
        Uy_fit = beta(5) + (beta(6) .* depth) + (beta(7).*  depth.^2) + (beta(8) .* depth.^3); 
        
    case 10    
        Ux_fit = beta(1) + (beta(2) .* depth) + (beta(3) .* depth.^2) + (beta(4) .* depth.^3) + (beta(5) .* depth.^4);
        Uy_fit = beta(6) + (beta(7) .* depth) + (beta(8).*  depth.^2) + (beta(9) .* depth.^3) + (beta(10) .* depth.^4);  
        
    case 12    
        Ux_fit = beta(1) + (beta(2) .* depth) + (beta(3) .* depth.^2) + (beta(4)  .* depth.^3) + (beta(5)  .* depth.^4) + (beta(6)  .* depth.^5);
        Uy_fit = beta(7) + (beta(8) .* depth) + (beta(9).*  depth.^2) + (beta(10) .* depth.^3) + (beta(11) .* depth.^4) + (beta(12) .* depth.^5) ;     
        
    otherwise
        error('wrong beta length')
end

end