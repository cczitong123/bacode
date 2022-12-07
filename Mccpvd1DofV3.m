classdef Mccpvd1DofV3 < Mccpvd1Dof
    %override dynamics: the control u = [u3], u1 u2 is predefined
    
    properties
        ffu1 %feedforward u1 motor trajectory of M1
        ffu2 %feedforward u2 motor trajectory of M2
    end
    
    methods
        function obj = Mccpvd1DofV3()
            obj = obj@Mccpvd1Dof();
        end
        function xdot = dynamics(obj, x, u3, n)
            u1 = obj.ffu1(n);
            u2 = obj.ffu2(n);
            qdd = obj.accel( x(1:2,:), [u1;u2;u3] );
            xdot = [x(2); qdd];
        end
    end
    
end

