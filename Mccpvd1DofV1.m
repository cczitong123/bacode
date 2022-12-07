classdef Mccpvd1DofV1 < Mccpvd1Dof
    %override dynamics: the control is [du1 du2 u3]'
    % state x: [q, qd, q1, q2]
    
    properties
        umax = [6.8 6.8 ]'
        umin = [-6.8 -6.8]'
    end
    
    methods
        function obj = Mccpvd1DofV1()
            obj = obj@Mccpvd1Dof();
        end
        
        function xdot = dynamics(obj, x, v)
            qdd = obj.accel(x(1:2,:), [x(3:4,:); v(3)]);
            xdot = [x(2); qdd; v(1); v(2)];
        end
        
        function [traj] = reachFastIlqr(obj, x0, qf)
            % traj opt using ilqr for fast reaching
            
            
            
        end
    end
    
end

