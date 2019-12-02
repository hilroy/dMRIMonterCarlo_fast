load('S0.mat');
load('Phi1.mat');
load('Phi2.mat');

% % Physical Constants
% D = 2; %mu_m^2/ms
% gamma = 2.675e5; % gyromagnetic ratio, unit: rad/ms/T

% % Sequence Parameters (ST PGSE)
% Delta = 37.7; %ms, separation of the two pulses 
% delta = 31.7; %ms, duration of the gradient
% Ttot = delta + Delta; %ms, total time for diffusion encoding
% TGap = Delta - delta;
% g = 4e-8; %T/mu_m
% g_hat = [1;0;0]; % gradient direction
% 
% Rcyl = 4;%confining domain is a cylinder of radius 4 mu_m
% L = 1000; % fibre length, approximately infinite
% bdlayer = Rcyl*1e-4;
% Rrels = Rcyl*0.1; % reflection radius
% 
% % distance between the particle and the boundary of the cylinder
% dist = @(r) Rcyl - norm(r(1:2));
% % interior normal unit vector for r at the boundary
% nvec = @(r) [-r(1:2)/norm(r(1:2));0];

T = delta;
Nwalker = 100;
Ph_sample = zeros(3,Nwalker);

tic
for k = 1:Nwalker
    % initial position inside the cylinder
    U = rand(3,1);
    r = [Rcyl*sqrt(U(1))*cos(pi*U(2)); 
        Rcyl*sqrt(U(1))*cos(pi*U(2));
        L*(2*U(3)-1)];
    d = dist(r);
    phi = [0;0;0];
    telps = 0;
    steps = 0;
    collision = 0;
    while telps < T
        tau = gen_tau;
        u = RandDir;
        if d > bdlayer % Walker is not near the boundary
            tc = d^2/D;
            if tc*tau < T - telps % not the final jump
                dt = tc*tau;
                jump = d;
                
                dphi = d*tc*evalu(tau,Phi1)*u; %phase increment
            else % last jump
                dt = T - telps;
                t_norm = dt/tc;
                jump = d*gen_r(t_norm);
                
                dphi = d*tc*evalu(t_norm,Phi2)*randn*[1;1;1];
            end
        else % Walker is at the boundary
            tc = Rrels^2/D;
            dt = tc*tau;           
            
            n = nvec(r);
            cosangle = dot(u,n);
            tgt = u - cosangle*n;
            tgt = tgt/norm(tgt);
            
            dphi = tc*Rrels*(n*tau/4 + tgt*evalu(tau,Phi1));
            
            if cosangle < 0 % reflect direction
                u = u - 2*cosangle*n;
            end
            jump = Rrels;
            collision = collision + 1;
        end
        phi = phi + dt*r + dphi;
        
        r = r + jump*u; 
        telps = telps + dt;
        d = dist(r);
        steps = steps+1;
    end
    Ph_sample(:,k) = phi;
end

SimTime = toc




