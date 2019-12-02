function [r_new, t_new, ph_new] = FRWupdate_regular(r,t,ph,dist, T_max, isgradient)
    % Walker is NOT near the boundary
    load('Phi1.mat');
    load('Phi2.mat');
    load('PhysicalConstants.mat');
     % random dimensionless first exit time
    tau = gen_tau;                   
    u = RandDir;
    % time_unit
    t_unit = dist^2 / D;              
    
    % not the final jump
    if  t + tau * t_unit < T_max
        % time increment is random
        dt = tau * t_unit; 
        t_new = t + dt;
        % new position is on the inscribed sphere
        r_new = r + dist * u;
        
        if isgradient == 1
            % nonlinear component of the phase increment
            % = Phi1 * displacement * time unit
            dph_nl = evalu(tau,Phi1) * u * dist * t_unit; 
            ph_new = ph + dt * r + dph_nl;
        else
            ph_new = ph;
        end
    % final jump
    else 
        t_new = T_max;
        % dimensionless time increment
        tau = (T_max - t) / t_unit;
        % new position is a random point inside the inscribed sphere
        r_new = r + dist * gen_r(tau) * u; 
        
        if isgradient == 1
            % dphi_nl is a normal random variable with sigma proportional
            % to Phi2
            dph_nl = evalu(tau,Phi2)* randn * ones(3,1) * dist * t_unit;
            ph_new = ph + (T_max - t) * r + dph_nl;
        else
            ph_new = ph;
        end
    end
end

