classdef STsequence
    % Diffusion gradient :Stejskal Tanner Pulsed Gradient Spin Echo 
    properties
        % Sequence Parameters
        Delta; % unit = ms, separation of the two pulses 
        delta; % unit = ms, duration of the gradient
        Ttot;
        Tgap;
        strength; %T/mu_m; eg g = 4e-8
        direction; % gradient direction, unit vector
        bvalue; % without gamma
    end
    
    methods
        function dgradient = STsequence(Delta,delta,strength,direction)
            % constructor
            dgradient.Delta     = Delta;
            dgradient.delta     = delta;
            dgradient.Ttot      = Delta + delta;
            dgradient.Ttot      = Delta - delta;
            dgradient.strength  = strength;
            dgradient.direction = direction;
            bvalue = strength^2 * delta^2 * (Delta - delta/3);
        end
    end
end

