classdef CrossFibre
    % orthogonally-crossing fibres
    properties
        iscvx;    % is the geometry convex?
        d_lat;    % halved lattice side length
        axis;     % 4 axes(see schematics), represented by segment end points 
        R_out;    % outer radius
        VF_out    % volume fraction of the extra-axonal space 
        d_bd;     % boundary layer thickness
        R_rls;    % release radius
    end
    
    methods
        function lattice = CrossFibre(d, R_out)
            lattice.iscvx = false(1);
            lattice.d_lat = d;
            A = zeros(3,2,4);
            A(:,:,1) = d*[1 -1; 1 1; -1 -1];   % x1
            A(:,:,2) = d*[1 -1; -1 -1; -1 -1]; % x2
            A(:,:,3) = d*[-1 -1; 1 -1; 1 1];   % y1
            A(:,:,4) = d*[1 1; 1 -1; 1 1];     % y2
            lattice.axis   = A;
            lattice.R_out  = R_out;
            lattice.VF_out = 1 - (pi/4)*(R_out/d)^2; 
            lattice.d_bd   = R_out * 1e-4;
            lattice.R_rls  = R_out * 0.05;
        end
        
        function [d_out, n_out] = dist(lat, pt)
            % dist(pt, an extra-axonal point, boundary)
            % here pt are the local coordinates wrt to the lattice containing
            % the spin, if dist is 'small', return a nonzero normal vector
            ax = lat.axis;
            d  = lat.d_lat;  
            index = floor(((pt/d)+1)/2); % addressing periodicity
            pt = pt - 2*d*index; 
            sides = pt - reshape(ax,[3,8]);
            side1 = sides(:,1:2:end-1);
            side2 = sides(:,2:2:end);
            DistAll = vecnorm(cross(side1,side2))/(2*d) - lat.R_out;
            d_out = min(DistAll);
            n_out = zeros(3,1);   % default value
            if 0 <= d_out && d_out < lat.d_bd             
                closest = find(DistAll == d_out);
                axdir = diff(ax(:,:,closest),1,2);
                side  = pt - ax(:,1,closest);
                n_out = side - dot(side,axdir)* axdir/(4*d^2);
                n_out = n_out/norm(n_out);
            end
        end
        
        function [r_ini_out, N_out, N_in] = unif(lat, N_walker)
            % generate initial positions
            d = lat.d_lat;
            N_out = round(N_walker * lat.VF_out);
            r_ini_out = zeros(3,1,N_out);
            n_out = 0;
            while n_out < N_out  % A-R scheme.
                r = d * (2*rand(3,1)-1);
                [d_out, ~] = dist(lat,r);
                if d_out > 0
                    n_out = n_out + 1;
                    r_ini_out(:,1,n_out) = r;
                end
                % expect a lot of unwanted samples ...
            end
            N_in = N_walker - N_out;
        end
    end
end