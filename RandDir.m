function [ r_hat ] = RandDir()
%Generate a 3d random unit column vector
U = rand(2,1);
z = 2*U(1)-1;
x = sqrt(1-z^2)*cos(2*pi*U(2));
y = sqrt(1-z^2)*sin(2*pi*U(2));
r_hat = [x;y;z];
end
