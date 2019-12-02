function [ jump ] = gen_r(t)
% t is the normalised time spent in the last sphere 
% jump is the distance between the walker and the centre of the sphere 
load('R(t).mat');
if t < 0.02 
    jump = 2;
    sigma = sqrt(2*t);
    while jump > 1
        jump_vec = sigma*randn(3,1);
        jump = norm(jump_vec);
    end
elseif t > 0.2
   U = rand;
   % for U close to 0 or 1 the solutions are calculated from taylor expansion 
   if U < 0.001 
       x = (3*U/pi^2)^(1/3);
   elseif U > 0.999
       x = 1-sqrt(2*(1-U)/pi);
   else % Newton's method
       A = U*pi;
       x = pi/2;
       f = sin(x)-x*cos(x)-A;
       df = x*sin(x);
       err = abs(f);
       while err > 1e-7
           x = x - f/df;
           f = sin(x)-x*cos(x)-A;
           df = x*sin(x);
           err = abs(f);
       end
       jump = x/pi;
   end;
else
    jump = evalu(t,Rtable);
end
end

