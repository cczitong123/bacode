classdef Qbmove1MD
    %QBMOVE1D Summary of this class goes here
    %   Detailed explanation goes here
    % x: q, qd, q1, q2, q1d, q1d
    % u: 
    
    properties
        dimQ = 1 % dimension of joints
        dimX = 6 % dimension of states
        dimU = 3 % dimension of controls
        
        inertia = 0.03 %3.19e-3
        Df = 0.002 % vis friction
        
        umax = [ pi/3; pi/3; 1] ;
        umin = [-pi/3; -pi/3; 0] ;
        
        actuator
        
        Dmax = 0.00848 % max damping
        
        alpha_servo = 20
    end
    
    methods
        function obj = Qbmove1MD(varargin)
            obj.actuator = Qbmove();
        end
        
        function [ xdot ] = dynamics(obj, x, u)
            qddot = obj.cmptAcc(x, u);
            qdot = x(2,:);
            xdot1 = [qdot; qddot] ;

            % motor dynamics
            xdot2 = obj.motor_dynamics_2nd(x(3:end,:),u(1:2,:));
            
            xdot  = [xdot1;
                    xdot2];
        end
        
        function [xdot, xdot_x, xdot_u]= dynamics_with_jacobian_fd(obj, x, u)
            qddot = obj.cmptAcc(x,u);
            qdot = x(2,:);
            xdot1 = [qdot; qddot] ;

            % motor dynamics
            [xdot2, xdot2_x, xdot2_u] = obj.motor_dynamics_2nd(x(3:end),u(1:2));
            
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
                 xdot_x = [0, 1, zeros(1,4) ;
                           daccdx           ;
                           zeros(4,2), xdot2_x];
           
                 xdot_u = [zeros(1,3) ;
                           daccdu     ; 
                           xdot2_u,zeros(4,1)];
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
        
        function [ xdot, xdot_x, xdot_u ] = motor_dynamics_2nd(obj, x, u)
            p = obj.alpha_servo;
            
            A = [ 0, 0, 1, 0;
                0, 0, 0, 1;
                -p^2, 0, -2*p, 0;
                0, -p^2, 0, -2*p];
            B = [ 0, 0;
                0, 0;
                p^2, 0;
                0, p^2];
            xdot = A*x + B*u;
            
            %xdot = A*x + B*u - [0;0; tau_l1/model.actuator.J1; tau_l2/model.actuator.J2];
            
            if nargout > 1
                xdot_x  = A;  
                xdot_u  = B; 
            end
            
        end
    end
    
end

function J = get_jacobian_fd ( f, x )

delta=1e-6;
y = f(x);
J = zeros(length(y),length(x));
for i=1:length(x)
	dx = zeros(size(x)); dx(i) = delta;
    try 
        yp = f(x+dx);
    catch
        yp = NaN;
    end
    try
        ym = f(x-dx);
    catch
        ym = NaN;
    end
    if isnan(yp)
        J(:,i) = -ym/delta;
    elseif isnan(ym)
        J(:,i) = yp/delta;
    else
        J(:,i) = ((yp - ym)/(2*delta));
    end
    
end
end