function [ output ] = evalu( t, table )
% t is the normalized time, 0.001<t<10 
% table is an arrary of 2 columns 
% each row is a pair of tabulated values eg, (t,S0(t))
a = 1;
b = size(table,1);
p = floor((a+b)/2);
ta = table(a,1)-t;
tp = table(p,1)-t;
i = 1;
    while (p-a > 1) && (i < 1000)
        if ta*tp < 0
            b = p;
        else
            a = p;
        end
        p = floor((a+b)/2);
        tp = table(p,1)-t;
    end
    % a and p are the desired indices
    % table(t) is computed via linear interpolation of point 
    % (t(a),table(a)),(t(p),table(p)) 
    A = table(a,:);
    P = table(p,:);
    slope = (P(2)-A(2))/(P(1)-A(1));
    output = A(2)+slope*(t-A(1));
end

