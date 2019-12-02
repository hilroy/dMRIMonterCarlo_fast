function [ tau ] = gen_tau( )
load('Ftau.mat');
% generate a random, dimensionless first-exit time, tau
U = rand();
FtauSmall = 1.45141496610268e-06;
FtauLarge = 0.999975026124371;

% when U < FtauSmall solve U = 2*exp(-1/(4*t))/sqrt(pi*t) using Netwon's
% method. The transendental equation can be written as
% f(t) = t*ln(t)+A*t+1/2, where A = 2*ln(U*sqrt(pi)/2) 
if U < FtauSmall
    A = 2*log(U*sqrt(pi)/2);
    t = 1; % initial guess;
    f = t*log(t)+A*t+1/2;
    err = abs(f);
    df = log(t)+1+A;
    while err > eps
        t = t - f/df;
        f = t*log(t)+A*t+1/2;
        err = abs(f);
        df = log(t)+1+A;
    end
    tau = t;
elseif U < FtauLarge 
    % for FtauSmall <= U < FtauLarge, we look up the table for Ftau. Using
    % a binary search scheme we look for the two consecutive row numbers 
    % in newFtau s.t. their newFtau values sandwich U
    a = 1;
    b = 2500;
    p = floor((a+b)/2);
    fa = newFtau(a,2)-U;
    fp = newFtau(p,2)-U;
    i = 1;
    while (p-a > 1) && (i < 1000)
        if fa*fp < 0
            b = p;
        else
            a = p;
        end
        p = floor((a+b)/2);
        fp = newFtau(p,2)-U;
    end
    % a and p are the desired indices
    % inverse_Ftau(U) is computed via linear interpolation of point 
    % (t(a),Ftau(a)),(t(p),Ftau(p)) 
    A = newFtau(a,:);
    P = newFtau(p,:);
    slope = (P(2)-A(2))/(P(1)-A(1));
    tau = P(1) - (P(2)-U)/slope;
else % when U > FtauLarge, use approximation 1-U = 2*exp(-pi^2*t)
    tau = log(2/(1-U))/(pi^2);
end
end

