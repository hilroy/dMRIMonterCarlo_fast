classdef SingleFibre
    properties
        iscvx; 
        R_in;
        d_bd
        R_rls;
    end
    
    methods
        function fibre = SingleFibre(R)
            fibre.iscvx = true(1);
            fibre.R_in = R;
            fibre.d_bd = R*1e-4;
            fibre.R_rls = R/20;
        end
        
        function [d_in, n_in] = dist(fib, pt)
            pt_xy = pt(1:2);
            d_in = fib.R_in - norm(pt_xy);
            n_in = zeros(3,1);
            if d_in < fib.d_bd
                n_in = [- pt_xy/norm(pt_xy);0];
            end
        end
        
        function r_ini_in = unif(fib,N_in)
            % generate uniformly random initial positions inside cylinder
            U = rand(1, N_in);
            rho = fib.R_in * sqrt(U);
            U = rand(1, N_in);
            theta = exp(2*1i*pi*U);
            dir = [real(theta); imag(theta)];
            r_ini_in = bsxfun(@times,dir,rho);
            r_ini_in = reshape([r_ini_in; zeros(1,N_in)],[3,1,N_in]);
        end
    end
end

