classdef Qbmove1D
    %QBMOVE1D velocity controlled servomotors
    %   x: q, qd, q1, q2
    %   u: v1, v2, d
    
    properties
        actuator
        Dmax = 0.0848
        
        dimQ = 1 % dimension of joints
        dimX = 4 % dimension of states
        dimU = 3 % dimension of controls
        
        inertia = 0.03 %3.19e-3
        Df = 0.002 % vis friction
    end
    
    methods
        function obj = Qbmove1D(varargin)
            %QBMOVE1D Construct an instance of this class
            %   Detailed explanation goes here
            obj.actuator = Qbmove();
        end
        
        function [ xdot ] = dynamics(obj, x, u)
            qddot = obj.cmptAcc(x, u);
            qdot = x(2,:);
            xdot1 = [qdot; qddot] ;

            % motor dynamics
            xdot2 = [u(1,:); u(2,:)];
            
            xdot  = [xdot1;
                    xdot2];
        end
        
        function [xdot, xdot_x, xdot_u]= dynamics_with_jacobian_fd(obj, x, u)
            qddot = obj.cmptAcc(x,u);
            qdot = x(2,:);
            xdot1 = [qdot; qddot] ;

            % motor dynamics
            %[xdot2, xdot2_x, xdot2_u] = obj.motor_dynamics_2nd(x(3:end),u(1:2));
            xdot2 = [u(1,:); u(2,:)];
            xdot  = [xdot1;
                     xdot2 ];
            % current motor positions
            %m     = [x(3);x(4)] ; 
            if nargout > 1
            % Compute xdot_x, xdot_u using finite differences
                %delta = 1e-6;
                %dimX = size(x,1);
                %dimU = size(u,1);
    
            % Compute derivative of acceleration w.r.t. x
                 %daccdx = zeros(1,model.dimX);
                 f = @(x)obj.cmptAcc(x,u);
                 daccdx = get_jacobian_fd(f,x);
                 if iscolumn(daccdx), daccdx = daccdx'; end
                 % Compute derivative of acceleration w.r.t. u
                 %daccdu = zeros(1,model.dimU);
                 f = @(u)obj.cmptAcc(x,u);
                 daccdu = get_jacobian_fd(f,u);
                 if iscolumn(daccdu), daccdu = daccdu'; end
                 xdot_x = [0, 1, 0, 0 ;
                           daccdx           ;
                           0, 0, 0, 0;
                           0, 0, 0, 0];
           
                 xdot_u = [zeros(1,3) ;
                           daccdu     ; 
                           1, 0, 0;
                           0, 1, 0];
             end
        end
        
        function acc = cmptAcc(obj, x, u)
            acc = obj.torque_total(x, u)./obj.inertia;
        end
        
        function tau = torque_total(obj, x, u)
            tau = obj.actuator.torque(x(1,:), x(3,:), x(4,:) ) ...
                - obj.Df*x(2,:) - obj.Dmax*u(3,:);
        end
        
        function s = stiffness(obj, x)
            s = obj.actuator.stiffness(x(1,:), x(3,:), x(4,:));
        end
        
        function d = damping(obj, u)
            d = obj.Dmax*u(3,:);
        end
    end
end

