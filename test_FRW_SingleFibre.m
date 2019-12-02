axis = [0;0;1];
radius = 3;
fib = SingleFibre(axis, radius);

Nwalker = 1;
r_ini = unif(fib,Nwalker);

isgradient = 0;
T_max = 20;


PhaseSample = zeros(3,Nwalker);

tic
for n = 1 : Nwalker
    [~, dephase] = FRW(fib, r_ini(:,n), T_max, isgradient);
    PhaseSample(:,n) = dephase;
end
toc