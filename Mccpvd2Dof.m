classdef Mccpvd2Dof < Arm2Dof
    %MCCPVD2D 
    % x: q1 q2 qd1 qd2 
    % u: (theta1 theta2 u3)_1, (theta1 theta2 u3)_2
    properties (Access = private)
      name = 'mccpvd2'
    end
    properties
        %name = 'mccpvd2'
        actuator1
        actuator2
        
        dimQ = 2
        dimU = 6
        
        L = [0.3; 0.3]; 
        I = [0.0045; 0.0045]; 
        M = [0.15; 0.15]; %[1.59; 1.44];
        Lg = [0.15; 0.15];
        g = 0;
        
        viscous_friction = 0.01;
        coulomb_friction = 0;
        
    end
    
    methods
        function obj = Mccpvd2Dof()
            obj = obj@Arm2Dof();
            obj.actuator1 = ActMccpvd();
            obj.actuator2 = ActMccpvd();
        end
        function torque = tau(obj, q, qdot, u)
            tau1 = obj.actuator1.torque(q(1,:), qdot(1,:), u(1,:), u(2,:), u(3,:));
            tau2 = obj.actuator2.torque(q(2,:), qdot(2,:), u(4,:), u(5,:), u(6,:));
            torque = [tau1; tau2];
        end
        
        function qddot = qddot(obj, q, qdot, u)
            % u: (m1 m2 u3)_1 , (m1 m2 u3)_2
            tau = obj.tau(q, qdot, u);
            qddot = Arm2Dof.compute_qddot(q, qdot, tau, obj);
        end
        
        % x: q1 q2 qd1 qd2
        % u: u
        function xdot = dynamics(model, x, u)
            qddot = model.qddot( x(1:2,:), x(3:4,:), u);
            xdot = [x(3:4); qddot];
        end
        
        function x = endpoint(obj, q)
            x = Arm2Dof.endpoint(q, obj.L); 
        end
        
        function J = jacobian(obj, q)
            J = Arm2Dof.jacobian(q, obj.L);
        end
        
    end
    
end

