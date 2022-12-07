classdef Mccpvd1DofV2 < Mccpvd1Dof
    %override dynamics: the control u = [u2 u3]', u1 is predefined
    
    properties
        ffu1
        
    end
    
    methods
        function obj = Mccpvd1DofV2()
            obj = obj@Mccpvd1Dof();
        end
        function xdot = dynamics(obj, x, u23, n)
            u1 = obj.ffu1(n);
            qdd = obj.accel( x(1:2,:), [u1; u23] );
            xdot = [ x(2); qdd];
        end
    end
    
end

