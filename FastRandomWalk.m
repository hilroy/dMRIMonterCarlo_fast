function spin_fin = FastRandomWalk(geometry,spin,N_walker,t_max,isgrad,const)
   % Fast Random Walk algorithm
   spin_fin = spin;
   for n_w = 1 : N_walker
       t_rem = t_max;
       r = spin(n_w).pos;
       ph = spin(n_w).phi;
       while t_rem > 0
           [d,nv] = dist(geometry,r);
           tau = gen_tau(const);
           u = RandDir;
           if d > geometry.d_bd % non-boundary case
               t_unit = d^2/const.D;
               dt = t_unit * tau;
               if  dt < t_rem
                   t_rem = t_rem - dt;     % r_new is on the inscribed sphere
                   r = r + d * u;
                   if isgrad ~= 0              % 'nl' stands for nonlinear
                       dph_nl = evalPhi1(const, tau)* d * t_unit * u;
                       ph = ph + (dt * r + dph_nl) * isgrad;
                   end
               else
                   dt = t_rem;
                   t_rem = 0;
                   tau = t_rem / t_unit;
                   r = r + d*gen_r(tau,const)*u; % r_new is inside the inscribed sphere
                   if isgrad ~= 0
                       % dphi_nl ~ N(0,sigma), where sigma is proportional to Phi2
                       dph_nl = evalPhi2(const, tau)* randn * d * t_unit * ones(3,1);
                       ph = ph + (dt * r + dph_nl) * isgrad;
                   end
               end
           else % boundary case
               R_rls = geometry.R_rls; 
               t_unit = R_rls^2/const.D;
               if geometry.iscvx == true(1)
                   while abs(dot(u,nv)) < 1/40 % avoid spins 'leaking' out!
                       u = RandDir;
                   end
               end
               % As R_rls << R_in, we assume the walker is always able to reach the
               % release hemisphere. 
               dt = tau * t_unit;
               t_rem = t_rem - dt;
               cosangle = dot(u,nv);
               tgt = u - cosangle * nv;
               tgt = tgt / norm(tgt);
               if isgrad ~= 0 
                   dph_nl = (nv*tau/4+tgt*evalPhi1(const, tau)) * R_rls * t_unit;
                   ph = ph + (dt * r + dph_nl) * isgrad;         
               end
               if cosangle < 0
                   u = u - 2*cosangle*nv;
               end
               r = r + u*R_rls;
           end
       end
       spin_fin(n_w).pos = r;
       spin_fin(n_w).phi = ph;
       if rem(n_w,10) == 0
           fprintf('%d out of %d iterations completed.\n', n_w, N_walker)
       end
   end
end
