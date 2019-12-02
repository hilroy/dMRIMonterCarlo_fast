classdef SingleFibre
    properties
        axis; % unit vector, default z+ direction
        radius;
        length; % length >> radius
        bdylayer;
        release;
    end
    
    methods
        function fibre = SingleFibre(axis, radius)
            % constructor
            if (nargin > 0) 
            fibre.axis = axis;
            fibre.radius = radius;
            fibre.bdylayer = radius*1e-4;
            fibre.length = 1000*radius;
            fibre.release = radius/20;
            end
        end
        
        function dist = distance(fibre, point)
            % if point is inside the fibre
            % calculate the distance between point and the cylinder
            a = fibre.axis;
            proj_axis = dot(point, a);
            dist = fibre.radius - norm(point - proj_axis*a) ;
            if dist < 0
                disp("point is outside the cylinder!");
            end
        end
        
        function nvector = normal(fibre, point)
            % if point is within the boundary layer
            % calculate the interior unit normal vector at point
            nvector = zeros(3,1);
            d = distance(fibre, point);
            if  d < fibre.bdylayer && d >= 0
                a = fibre.axis;
                n = point - dot(point, a)*a;
                nvector = -n/norm(n);
            else
                disp("point is far away from the boundary!");
            end
        end
        
        function ini_pos = unif(fibre,Nwalker)
            % generate uniformly random initial positions inside the
            % cylinder
            r = fibre.radius;
            ini_pos = zeros(3,Nwalker);
            U = rand(3,Nwalker);
            ini_pos(1:2,:) = r*sqrt(U(1,:)).*[cos(2*pi*U(2,:));sin(2*pi*U(2,:))];
            ini_pos(3,:) = fibre.length*(U(3,:)-1/2);
        end
        
        function [r_fin, dephase] = FRW(fibre, r_ini, T_max, isgradient)
            % fast random walk algorithm for 1 walker
            thic = fibre.bdylayer;
            R_rls = fibre.release;
            
            t_elapsed = 0;
            r = r_ini;
            dephase = zeros(3,1);
            
            while t_elapsed < T_max
                dist = distance(fibre, r);
                if dist > thic
                    [r, t_elapsed, dephase] = FRWupdate_regular(r,t_elapsed,dephase,dist, T_max, isgradient);
                else
                    n = normal(fibre, r);
                    [r, t_elapsed, dephase] = FRWupdate_boundary(r, t_elapsed, dephase, n, R_rls, isgradient);
                end
            end
            r_fin = r;
        end
    end
end

