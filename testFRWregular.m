r = zeros(3,1);
t = 0;
ph = zeros(3,1);

dist = 1.5;
T_max = 0.2;
isgradient = 1;

[r_new, t_new, ph_new] = FRWupdate_regular(r,t,ph,dist, T_max, isgradient);