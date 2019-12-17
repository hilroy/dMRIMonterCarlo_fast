classdef Tables_FRW
    % tabulated values of F_tau, Phi_1, Phi_2, D and gamma
    properties
        D;
        gam;
        F_tau;
        Phi_1;
        Phi_2;
        R_rms;
    end
    methods
        function const = Tables_FRW()
            load('PhysicalConstants.mat','D','gamma')
            load('Ftau.mat','newFtau')
            load('Phi1.mat','Phi1');
            load('Phi2.mat','Phi2');
            load('R(t).mat','Rtable');
            const.D = D;
            const.gam = gamma;
            const.F_tau = newFtau(2:135,:);
            const.Phi_1 = Phi1(1:800,:);
            const.Phi_2 = Phi2;
            const.R_rms = Rtable(1:200,:);
        end
        
        function tau = gen_tau(const)
            U = rand;
            while U == 1
                U = rand;
            end
            if U < 1e-5
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
            elseif U < 0.99
                tau = interp1(const.F_tau(:,2),const.F_tau(:,1),U);
            else
                tau = log(2/(1-U))/(pi^2);
            end
        end
        
        function r_rms = gen_r(tau, const)
            if tau < 0.02 
                r_rms = 2;
                sigma = sqrt(2*tau);
                while r_rms > 1
                    dr = sigma*randn(3,1);
                    r_rms = norm(dr);
                end
            elseif tau < 0.2
                r_rms = interp1(const.R_rms(:,1),const.R_rms(:,2),tau);
            else
                U = rand;
                if U < 0.001 
                    r_rms = (3*U/pi^2)^(1/3); % taylor expansion 
                elseif U > 0.999
                    r_rms = 1-sqrt(2*(1-U)/pi);
                else % Newton's method
                    A = U*pi;
                    x = pi/2;
                    f = sin(x)-x*cos(x)-A;
                    df = x*sin(x);
                    err = abs(f);
                    while err > 1e-7
                        x = x - f/df;
                        f = sin(x)-x*cos(x)-A;
                        df = x*sin(x);
                        err = abs(f);
                    end
                    r_rms = x/pi;
                end
            end
        end
        
        function output = evalPhi1(const, tau)
            if tau < 0.05
                output = tau*(1-2*tau)/2;
            elseif tau < 0.8
                output = interp1(const.Phi_1(:,1),const.Phi_1(:,2),tau);
            else
                output = 0.75/pi^2;
            end
        end
           
        function output = evalPhi2(const, tau)
            if tau < 1e-3
                output = sqrt(2*tau^3/3);
            else
                output = interp1(const.Phi_2(:,1),const.Phi_2(:,2),tau);
            end
        end
    end
end


