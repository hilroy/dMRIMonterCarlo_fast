r = zeros(3,1);
n = [0;0;1];
ph = zeros(3,1);
t = 0;
R_rls = 1.5;
isgradient = 1;
FRWupdate_boundary(r,t,ph,n, R_rls, isgradient)

[r_new, t_new, ph_new] = FRWupdate_boundary(r,t,ph,n, R_rls, isgradient);