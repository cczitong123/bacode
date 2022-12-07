classdef ActMccpvd
    %ACTMCCPVD Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %%%% From servomotor datasheet
        %%% no load speed: 6.16 at 6V, 7.48 at 7.4V
        %%% no load current: 0.36,     0.44
        %%% stall torque: 2.8439285, 3.4323275
        %%%% Given the datasheet
        % we want to constraint the current < 3A, 
        % max torque:     1.14
        % max accel:      115
        %%%% Geometry parameters
        B = 0.036
        C = 0.135
        A0 = 0.099 % C - B
        Ks = 394 %231 % spring constant
        r = 0.015
        gear = 1; %3/4 % gear ratio between servo1 and link axis
        %damping_range = [0, 0.00848]; % 0.00848
        damping_range %= [0, 0.00848]; % 0.00848
        
        %%%% servos parameters
        J1 = 0.0099;
        K1 = 0.3811; % torque constant
        D1 = 0.023; %0.2009;
        R1 = 0.822;
        J2 = 0.0099;
        K2 = 0.3811;
        D2 = 0.023; %0.2009;
        R2 = 0.822;
        
        %b = 0.0007382;
        %b = 0;
        b = 0.0002;
        
        %%%% - damping motor - %%%%
        motor_inertia = 4.6*10^(-6);
        Id = 0.00184
        gear_d = 40;
        Kd = 0.0212; %torque constant
        Rd = 21.2;
        Rl = 25.3;
        
        ratio_load = 21.2/25;
        max_damping_db
        max_damping
        max_rege_damping
        
        u_max_regedamp
    end
    
    methods
        function obj = ActMccpvd(varargin)
            if nargin > 0
                param = varargin{1};
                if isfield(param,'ratio_load'), obj.ratio_load = param.ratio_load; end
                if isfield(param,'gear_d'), obj.gear_d = param.gear_d ; end
                if isfield(param,'Kd'), obj.Kd = param.Kd ; end
                if isfield(param,'J1'), obj.J1 = param.J1; end
                if isfield(param,'J2'), obj.J2 = param.J2; end
                if isfield(param,'K1'), obj.K1 = param.K1; end
                if isfield(param,'K2'), obj.K2 = param.K2; end
                if isfield(param,'R1'), obj.R1 = param.R1; end
                if isfield(param,'R2'), obj.R2 = param.R2; end
                if isfield(param,'Ks'), obj.Ks = param.Ks; end
                if isfield(param,'Rl')
                    obj.Rl = param.Rl; obj.ratio_load = obj.Rl/obj.Rd; 
                end
            end
            
            
            
            
            obj.max_damping = obj.Kd^2 * obj.gear_d^2/obj.Rd ;
            obj.max_rege_damping = obj.Kd^2 * obj.gear_d^2/(obj.Rd+obj.Rl);
            obj.u_max_regedamp = 1/(1 + obj.ratio_load);
            %obj.u_max_regedamp = 0.5;
            obj.Id = obj.motor_inertia*obj.gear_d^2;
        end
        
        %variable damping
        function d = damping(obj, u)
             
             ratio = obj.ratio_load;
             
             DC1 = max(min(u/obj.u_max_regedamp, 1),0);
             DC2 = max(min(( u - obj.u_max_regedamp )./( 1 - obj.u_max_regedamp ), 1),0);
             d = obj.max_rege_damping*DC1 + obj.max_rege_damping*DC2*ratio;
             
            %d = obj.max_damping*u3;
        end
        
       
        
        function power = power_rege(obj, qdot, u)
            ratio = obj.ratio_load;
            alpha = ratio/(1+ratio);
            
            DC1 = max(min(u/obj.u_max_regedamp, 1),0);
            DC2 = max(min(( u - obj.u_max_regedamp )./( 1 - obj.u_max_regedamp ), 1),0);
            power = obj.max_rege_damping*alpha*(qdot.^2).*(DC1 - DC2);
        end
        
        function p = p_damp_inputelec(obj, qdot, u)
            p = obj.gear_d^2*obj.Kd^2*qdot^2*obj.transm(u)^2/(obj.Rd+obj.Rl);
        end
        function p = p_damp_inputmech(obj, qdot, u)
            p = obj.gear_d^2*obj.Kd^2*qdot^2*obj.transm(u)/(obj.Rd+obj.Rl);
        end
        function torque_damp = torque_damping(obj,qdot, dc)
            % damping torque
            % note size(x,2) == size(u,2)
            torque_damp = qdot.*obj.damping(dc);
        end
        
        
        
        function torque = torque(obj, q, qdot, m1, m2, dc)
            torque = obj.torque_spring(q, m1, m2) - obj.torque_damping(qdot,dc);
        end
        
        function torque = torque_spring(obj, q, m1, m2)
            %sprdis = model.spring_displacement(x,u);
            %spring_force = sprdis*model.spring_constant;
            phi = m1/obj.gear - q;
            m2 = max(0,m2);
            A = sqrt( obj.C^2 + obj.B^2 - 2*(obj.B*obj.C).*cos( phi ));
            %m2 = max( (obj.A0-A)/obj.r ,0);
            torque = obj.Ks*(obj.B*obj.C).*sin(phi).*...
                (1+ (obj.r*m2 - obj.A0)./A );
        end
        function tau = torque_load1(obj, q, theta1, theta2)
            % tau = - dU/dtheta1
            tau = -obj.torque_spring(q, theta1, theta2)/obj.gear;
        end
        function tau = torque_load2(obj, q, theta1, theta2)
            % tau = - dU/dtheta2
            phi = theta1/obj.gear - q;
            A = sqrt( obj.C^2 + obj.B^2 - 2*(obj.B*obj.C).*cos( phi ));
            tau = -obj.Ks*(A + obj.r*theta2 - obj.A0)*obj.r;
        end
        function force = spring_force(obj, q, theta1, theta2)
            sprlength = obj.spring_displacement(q, theta1, theta2);
            force = sprlength*obj.Ks;
        end
        function sprlength = spring_displacement(obj, q, theta1, theta2)
            
            A = sqrt(obj.C^2 + obj.B^2 - 2*obj.B*obj.C*cos(theta1/obj.gear - q) );
            sprlength = A + obj.r*theta2 - obj.A0; 
        end
        
        function k = stiffness(obj, q, theta1, theta2)
            
            A = sqrt( obj.C^2 + obj.B^2 - 2*obj.B*obj.C.*cos( theta1/obj.gear - q ) );
            phi = theta1/obj.gear - q;
            k=  obj.Ks*obj.B*obj.C*cos(phi).*( 1 + (obj.r*theta2 - obj.A0)./A ) - ...
                obj.Ks*(obj.B*obj.C*sin(phi)).^2.*( obj.r*theta2 - obj.A0)./A.^1.5 ;
        end
        
        
    end
    
    methods (Static)
        function tr = transm(u)
            % u :- u3
            tr = u;
        end
        function mddot = motor_accel(m, mdot, u)
            beta = 25;
            mddot = beta^2*(u - m) - 2*beta*mdot;
        end
    end
    
end

