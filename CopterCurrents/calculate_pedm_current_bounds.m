function [U_vertices, Z_vertices, Umin, Umax] = calculate_pedm_current_bounds(PEDM,eps_PEDM,verbose,zvect,percentVal)

%PEDM - either a function handle or a set of polynomial coffiecients
%defining the PEDM profile as a function of depth.

if isa(PEDM,'function_handle')
    U_PEDM_z = PEDM(zvect);
else
    U_PEDM_z = polyval(pPEDM,zvect);
end

Umin = U_PEDM_z;
Umax = U_PEDM_z;

for i = 1:numel(verbose)
   if verbose(i).eps_PEDM < eps_PEDM*(1+percentVal)
       U_iz = polyval(verbose(i).pPEDM,zvect);
       
       indsMax = U_iz > Umax;
       indsMin = U_iz < Umin;
       
       Umax(indsMax) = U_iz(indsMax);
       Umin(indsMin) = U_iz(indsMin);  
      

       
   end
        
end


Z_vertices = [zvect,zvect(end:-1:1)];
U_vertices = [Umin,fliplr(Umax)];


end