function [r_new, t_new, ph_new] = FRWupdate_boundary(r,t,ph,n, R_rls, isgradient)
     % Walker is near the boundary
     load('Phi1.mat');
     load('PhysicalConstants.mat');
     
     tau = gen_tau;
     u = RandDir;
     t_unit = R_rls^2 / D;
     
     % As R_rls is small compared to the geometry's characteristic length,
     % we do not expect updated time to exceed T_max for the boundary jump
     dt = tau * t_unit;
     t_new = t + dt;
     
     cosangle = dot(u,n);
     tgt = u - cosangle * n;
     tgt = tgt / norm(tgt);
     
     if isgradient == 1       
         dph_nl = (n * tau / 4 + tgt * evalu(tau,Phi1)) * R_rls * t_unit;
         ph_new = ph + dt * r + dph_nl;
     else
         ph_new = ph;
     end
     
     if cosangle < 0 
         % reflect direction
         u = u - 2*cosangle*n;
     end
     r_new = r + u * R_rls;
end

