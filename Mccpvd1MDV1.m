classdef Mccpvd1MDV1 < Mccpvd1MD
    %MCCPVD1MDV1 control with predefined u1 and q_m1
    % x: [ q, qdot, q_s, qdot-s]
    
    properties
        % F: Tau_a
        df1n % dq
        df2n % dqdot
        df3n % dqm1
        df4n % dqm2
        df5n % du3
    end
    
    methods
        function obj = Mccpvd1MDV1()
            obj = obj@Mccpvd1MD();
            obj.dimQ = 1; % dimension of joints
            obj.dimX = 4; % dimension of states
            obj.dimU = 2; % dimension of controls
            
            B = obj.actuator.B;
            C = obj.actuator.C;
            A0 = obj.actuator.A0; % C - B
            Ks = obj.actuator.Ks; %231 % spring constant
            r = obj.actuator.r;
            Dm = obj.actuator.max_damping;
            
            % q, qdot, th1, th2, th3
            syms F(x1,x2,x3,x4,x5)
            
            F(x1,x2,x3,x4,x5) = Ks*B*C*sin(x3-x1)*(1+ (r*x4 - A0)/sqrt(B^2+C^2-2*B*C*cos(x3-x1)) ) - Dm*x5*x2;
            
            
            df1 = diff(F,x1);
            df2 = diff(F,x2);
            df3 = diff(F,x3);
            df4 = diff(F,x4);
            df5 = diff(F,x5);
            
            obj.df1n = matlabFunction(df1);
            obj.df2n = matlabFunction(df2);
            obj.df3n = matlabFunction(df3);
            obj.df4n = matlabFunction(df4);
            obj.df5n = matlabFunction(df5);
        end
        function [xdot, xdot_x, xdot_u] = dynamics(obj, x, u23, t, u1ref, tref)
            idx = zeros(size(t));
            for i=1:length(idx)
                idx(i) = find(t(i)>=tref, 1, 'last');
            end
            u1 = u1ref(idx);
            %u = [u1; u23];
            
            dx1 = x(2, :);
            dx3 = x(4, :);
            alpha = obj.alpha_servo;
            dx4 = alpha^2*( u23(1,:) - x(3,:) ) - 2*alpha*x(4,:);
            
            dx2 = obj.cmptAcc(x, u23, u1);
            
            xdot  = [ dx1; dx2; dx3; dx4 ];
            M = obj.inertia;
            fv = obj.Df;
            if nargout > 1
                xd_x_1 = [0 1 0 0];
                xd_x_2 = [obj.df1n(x(1,:), x(2,:), u1, x(3,:), u23(2,:))/M, ...
                         (obj.df2n( x(1,:), x(2,:), u1, x(3,:), u23(2,:))- fv)/M,...
                         (obj.df4n(x(1,:), x(2,:), u1, x(3,:), u23(2,:)))/M,...
                         0];
                xd_x_3 = [0 0 0 1];
                xd_x_4 = [0 0 -alpha^2 2*alpha];
                xdot_x = [xd_x_1; xd_x_2; xd_x_3; xd_x_4];
                
                xdot_u = [0, 0;...
                         0, obj.df5n( x(1,:), x(2,:), u1, x(3,:), u23(2,:) )/M;...
                         0, 0;...
                         alpha^2, 0];
            end
        end
        function acc = cmptAcc(obj, x, u, u1ref)
           acc = obj.torque_total(x, u, u1ref)./obj.inertia;
        end
        function torque_total = torque_total(obj, x, u, u1ref)
            torque_total = obj.actuator.torque(x(1,:), x(2,:), u1ref, x(3,:), u(2,:))...
                - obj.Df.*x(2,:) - sin(x(1,:))*obj.mass*obj.gravity_constant*obj.com;
        end
    end
    
end

