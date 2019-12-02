classdef CrossFib3DCell
    % define the periodic structure for orthogonally-crossing fibres
    % only the space BETWEEN the fibres are considered in this class
    properties
        length; % half side length of the cube
        radius; % radius of cylindrical fibres
        bdylayer; % boundary layer thickness
        axis;       
    end
    
    methods
        function cell = CrossFib3DCell(length, radius)
            % class constructor
            if (nargin > 0)
            cell.length = length;
            cell.radius = radius;
            cell.bdylayer = radius*1e-4;
            A = zeros(3,2,4);
            A(:,:,1) = length*[1 -1; 1 1; -1 -1];
            A(:,:,2) = length*[1 -1; -1 -1; -1 -1];
            A(:,:,3) = length*[-1 -1; 1 -1; 1 1];
            A(:,:,4) = length*[1 1; 1 -1; 1 1];
            cell.axis = A;
            end
        end
        
        function [dist, closest] = distance(cell, point)
            % calculate the distance between cell and point
            % also, return the closest axis
            A = cell.axis;
            R = cell.radius;
            DistAllAx = zeros(1,4);
            for k = 1:4
                DistAllAx(k) = DistP2L(A(:,:,k),point) - R;
            end
            dist = min(DistAllAx);
            closest = find(DistAllAx == dist);         
        end
        
        function nvector = normal(cell, point)
            % if point is within the boundary layer of the cell
            % calculate the normal vector of the surface at point 
            % which is outward-pointing and has length = 1.       
            
            [dist, closest] = distance(cell, point);
            nvector = zeros(3,1);
            A = cell.axis(:,:,closest);
            
            if dist < cell.bdylayer
                a = point - A(:,1);
                b = A(:,2) - A(:,1);
                aprojb = (dot(a,b)/dot(b,b))*b;
                aperpb = a - aprojb;
                nvector = aperpb/norm(aperpb);
            end        
        end
        
        function ini_pos = unif(cell, Nwalker)
            % generate number = Nwalker of initial positions
            ini_pos = zeros(3,Nwalker);
            l = cell.length;
            n = 1;
            while n <= Nwalker
                r = l*(2*rand(3,1)-1);
                d = distance(cell,r);
                if d > 0
                    ini_pos(:,n) = r;
                    n = n + 1;
                end
            end
        end
        
        function r_tran = translate(cell, r)
            % r : position, 3 by 1
            % if r is outside of the cell, r_tran is the translated r
            % this is to address periodic boundary position
            l = cell.length;
            d = 2*l;
            if r(1) > l
                r(1) = r(1) - d;
            elseif r(1) < -l
                r(1) = r(1) + d;
            end
            if r(2) > l
                r(2) = r(2) - d;
            elseif r(1) < -l
                r(2) = r(2) + d;
            end
            if r(3) > l
                r(3) = r(3) - d;
            elseif r(1) < -l
                r(1) = r(3) + d;
            end
        end
        
        function [r_fin, ]
        
end

