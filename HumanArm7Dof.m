classdef HumanArm7Dof
    % HUMANARM7DOF 7 Dof 3 link human arm model
    %   Detailed explanation goes here
    
    properties
        
        %x % states vector
        
    end
    
    methods
        
        % joint torques
        function tau = torque(q, qdot, a)
            
            
            tau = M*T;
        end
        
        function M = moment_arm(q)
        
        end
        
        function T = muscle_tension(q, qdot, a)
        
        end
        
        function adot = muscle_dynamics(a, u)
            alpha = 40; % 40 msec
            adot = (u - a)/alpha;
        end
    
end

