classdef Mccpvd1MD < handle
    %   MACCEPA-VD
    %   Detailed explanation goes here
    %   1DoF MACCEPAVD
    %   x: q ; dq ; theta1 ; theta2; dtheta1; dtheta2;
    %   
    %   
    
    properties
        robot_type = 'Maccepavd1DOF-2nd Motor Dyn';
        
        cdt = 0.02 % controller time step
        
        dimQ = 1 % dimension of joints
        dimX = 6 % dimension of states
        dimU = 3 % dimension of controls
        
        
        %%%% physical design
        % physical joint constriants
        % note that the physical constriant of joint is the rechable point
        % not the controllable point (for EP motor)
        qmax =  135*pi/180; %
        qmin = -135*pi/180; %
        
        link_length = 0.15; % link length
        link_mass = 0.09;
        com = 0.085;
        
        servo1_mass = 0.09;
        
        mass = 0.18; % servo1_mass + link_mass
        
%         rege_ratio = 0;
        
        
        %%%% ---- physical design
        
        %%%% servo motor specification, reflected at gearbox output shaft
%         inertia_m1 = 0.0099;
%         K_m1 = 0.3811; % torque constant
%         D_m1 = 0.2009 ;
%         R_m1 = 0.8222;
%         inertia_m2 = 0.0099;
%         K_m2 = 0.3811;
%         D_m2 = 0.2009;
%         R_m2 = 0.8222;
        
        alpha_servo = 20; % Fan: fit from data
        %%%%
        
        %%%% dynamical properties
        % Fan: estimated by data  %calculated inertia: 0.00135
        % 
        inertia = 3.6e-3 %3.19e-3
        inertia_l = 0.00135
        % frictions
        % viscous friction
        % Fan: 0.0022 estimated by data;
        Df = 0.0077
        % coulomb friction
        coulomb_friction = 0;
        % gravity constant
        gravity_constant = 0;
        
        %%%%
        %%%% control input limits
        % u1 and u2 are defined in rad; round to int and below the physical
        % limit to protect the servo
        umax = [ pi/3; pi/2; 1] ;
        umin = [-pi/3; pi/6; 0] ;
        %%%%
        
        
        
        %%%% function handlers
        fnMotorDyn
        dynx
        dynu
        %%%%
        
        %%%% symbolic
        symDyn
        symDynx
        symDynu
        
        
        %%%%
        % Including another object and calling its function may be slower
        % in Matlab for simulating the dynamics. Consider define the actuator model
        % directly within this class.
        actuator
        %%%%
        
        
    end
    
    properties (Access = private)
        dynxfull
        dynufull
    end
    
    methods
        function model = Mccpvd1MD(varargin)
            
            % gear ratio from theta1 to q
%             gear1 = 3/4; %Nin = 44, Nout = 33; 
%         
%             spring_constant = 231;
%             lever_length = 0.036; %B
%             pin_displacement = 0.135; %C
%             
%             drum_radius = 0.015;
            %variable damping range
            % 110255 : 0.00848
%             damping_range = [0, 0.00848];
            if nargin >= 2
                model.actuator = ActMccpvd(varargin{2});
            else
                model.actuator = ActMccpvd();
            end
            if nargin > 0
                param = varargin{1};
                if isfield(param,'inertia_l')
                    model.inertia_l = param.inertia_l; 
                    model.inertia = model.inertia_l + model.actuator.Id;
                end
                if isfield(param,'inertia'), model.inertia = param.inertia; end
                if isfield(param,'Df'), model.Df = param.Df; end
                if isfield(param,'has_gravity'), model.gravity_constant = 9.8; end
                %if isfield(param,'gear'), obj.gear_d = param.gear_d ; end
            end
            
            %model.inertia = model.inertia_l + model.actuator.Id;
            %model.mass = model.link_mass + model.servo1_mass;
            
            %model = model.init_symfuns();

        end
        
%         function [model] = init_symfuns(model)
%             %%%% setup symbolic dynamics and compute jacobian functions
%             xv=sym('x',[6 1]);
%             uv=sym('u',[3 1]);
%             model.symDyn = symfun(model.dynamics(xv,uv),[xv;uv]);
%             model.symDynx = jacobian(model.symDyn,xv);
%             model.symDynu = jacobian(model.symDyn,uv);
%             model.dynxfull = matlabFunction(model.symDynx);
%             model.dynufull = matlabFunction(model.symDynu);
%             model.dynx =@(x,u)model.dynxfull (x(1),x(2),x(3),x(4),x(5),x(6),u(1),u(2),u(3));
%             model.dynu =@(x,u)model.dynufull(x(1),x(2),x(3),x(4),x(5),x(6),u(1),u(2),u(3));
%             %model.dynx = matlabFunction(model.dynx);
%             %model.dynu = matlabFunction(model.dynu);
%             %%%%
%             
%         end
        
        % state-space forward dynamics
        function [xdot]= dynamics(model, x, u)
            qddot = model.cmptAcc(x, u);
            qdot = x(2,:);
            xdot1 = [qdot; qddot] ;

            % motor dynamics
            xdot2 = model.motor_dynamics_2nd(x(3:end,:),u(1:2,:));
            
            xdot  = [xdot1;
                    xdot2];

        end
        
        function [xdot, xdot_x, xdot_u]= dynamics_with_jacobian_fd(model, x, u)
            qddot = model.cmptAcc(x,u);
            qdot = x(2,:);
            xdot1 = [qdot; qddot] ;

            % motor dynamics
            [xdot2, xdot2_x, xdot2_u] = model.motor_dynamics_2nd(x(3:end),u(1:2));
            
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
                f = @(x)model.cmptAcc(x,u);
                daccdx = get_jacobian_fd(f,x);
                if iscolumn(daccdx), daccdx = daccdx'; end
                % Compute derivative of acceleration w.r.t. u
                %daccdu = zeros(1,model.dimU);
                f = @(u)model.cmptAcc(x,u);
                daccdu = get_jacobian_fd(f,u);
                if iscolumn(daccdu), daccdu = daccdu'; end
                xdot_x = [0, 1, zeros(1,4) ;
                          daccdx           ;
                          zeros(4,2), xdot2_x];
                
                xdot_u = [zeros(1,3) ;
                          daccdu     ;
                          xdot2_u, zeros(4,1)];
             end
        end
        
        % incomplete
        % dynamics with analytical derivatives
        function [xdot, xdot_x, xdot_u] = dynamics_with_jacobian(model, x, u)
            xdot = model.dynamics(x,u);
            if nargout > 1
                xdot_x = model.dynx(x,u);
                xdot_u = model.dynu(x,u);
            end
        end
        
        function acc = cmptAcc(model,x,u)
           acc = model.torque_total(x,u)./model.inertia;
        end
        %%%% dynamics %%%%
        
        %%%%---- dynamical properties
        function k = stiffness(model, x)           
            k =  model.actuator.stiffness(x(1,:),x(3,:),x(4,:));
        end
        function dr = damping_ratio(model, x, u)
            b = model.Df + model.damping(u);
            k = model.stiffness(x);
            dr = b/(2*sqrt(k*model.inertia));
        end
        
        %variable damping
        function d = damping(model, u)
            % duty_circle: 0-1
            % linear on duty-circle
            d = model.actuator.damping(u(3,:));
        end
        %%%% ---- dynamical properties ---- %%%%
        
        %%%%---- torques ----%%%%
        % Spring Torque
        function torque = torque_spring(model, x,~)
            %sprdis = model.spring_displacement(x,u);
            %spring_force = sprdis*model.spring_constant;
            
            torque = model.actuator.torque_spring(x(1,:),x(3,:),x(4,:));
        end
        
        function force = spring_force(model, x, ~)
            
            %force = model.actuator.spring_force(x(1), x(3), x(4));
            B = model.actuator.B;
            C = model.actuator.C;
            A0 = model.actuator.A0;
            phi = x(3,:) - x(1,:);
            A = sqrt(B^2 + C^2 - 2*B*C*cos(phi));
            force = model.actuator.Ks*( A + model.actuator.r*x(4,:) - A0 );
        end
        % Variable Damping torque
        function torque_damp = torque_damping(model,x,u)
            % note size(x,2) == size(u,2)
            torque_damp = model.actuator.torque_damping(x(2),u(3));
        end
        % Actuator torque
        function torque = torque(model, x, u)
            % actuator torque
            torque = model.actuator.torque(x(1),x(2),x(3),x(4),u(3));
        end
        % Torque including friction, gravity
        function torque_total = torque_total(model, x, u)
            torque_total = model.actuator.torque(x(1,:),x(2,:),x(3,:),x(4,:),u(3,:))...
                - model.Df.*x(2,:) - sin(x(1,:))*model.mass*model.gravity_constant*model.com;
        end
        %%%%---- torques ----%%%%
        
        %%%%---- motor torques ----%%%%
        function tau_m = tau_m1(model, x, u)
            % tau_m of servo1 at servo's shaft
            % M ddtheta + D dtheta = tau_m + tau_l
            tau_l = model.actuator.torque_load1(x(1,:),x(3,:),x(4,:));
            p = model.alpha_servo;
            accel = p^2*(u(1,:)-x(3,:)) - 2*p*x(5,:);
            
            tau_m = - tau_l + model.actuator.J1*accel + model.actuator.D1*x(5,:) ;
            %tau_m = model.actuator.J1*accel + model.actuator.D1*x(5) ;
            %tau_m = tau_l;
        end
        
        function tau_m = tau_m2(model, x, u)
            % tau_m of servo2 at servo's shaft
            tau_l = model.actuator.torque_load2(x(1,:),x(3,:),x(4,:));
            p = model.alpha_servo;
            accel = p^2*(u(2,:)-x(4,:)) - 2*p*x(6,:);
            
            tau_m = - tau_l + model.actuator.D2*x(6,:) + model.actuator.J2*accel;
            %tau_m = model.actuator.D2*x(6) + model.actuator.J2*accel;
        end
        
        %%%% ---- motor torques
        
        %%%% power
        function out = sqrt_torques(model, x, u)
            t1 = model.actuator.torque_load1(x(1,:),x(3,:),x(4,:));
            t2 = model.actuator.torque_load2(x(1,:),x(3,:),x(4,:));
            out = t1.^2 + t2.^2;
        end
        function p = power_out(model, x, u)
            % to replace power_link
            tau = model.torque_spring(x,u);
            p = tau.*x(2,:);
        end
        function p = power_in(model, x, u)
            %tau_l1 = model.actuator.torque_load1(x(1,:),x(3,:),x(4,:));
            %tau_l2 = model.actuator.torque_load2(x(1,:),x(3,:),x(4,:));
            %p = max(-tau_l1.*x(5,:),0) + max(-tau_l2.*x(6,:),0);
            p1 = model.power_in1(x,u);
            p2 = model.power_in2(x,u);
            p1 = max(p1,0);
            p2 = max(p2,0);
            p = p1 + p2;
        end
        function p = power_in1(model, x, ~)
            tau_l1 = model.actuator.torque_load1(x(1,:),x(3,:),x(4,:));
            p = -tau_l1.*x(5,:);
        end
        function p = power_in2(model, x, ~)
            tau_l2 = model.actuator.torque_load2(x(1,:),x(3,:),x(4,:));
            p = -tau_l2.*x(6,:);
        end
        function p = power_damp(model,x, u)
            % energy dissipated via damping
            p = model.damping(u).*(x(2,:).^2 );
        end
        function p = power_rege(model, x, u)
            % power of regeneration
            p = model.actuator.power_rege(x(2,:),u(3,:));
        end
        function E = energy_spring(model, x, ~)
            l = model.actuator.spring_displacement(x(1,:),x(3,:),x(4,:));
            E = model.actuator.Ks*(l.^2)/2;
        end
        %%%% to-do: the following power calculation need reviewing
        function [power, p1, p2] = power_mech(model,x,u)
            tau_m1 = model.tau_m1(x,u);
            tau_m2 = model.tau_m2(x,u);
            p1 = tau_m1*x(5);
            p2 = tau_m2*x(6);
            power = p1 + p2;
            
            
        end
        
        function [power, p1_elec, p2_elec, p1_diss, p2_diss] = power_elec2(model, x, u)
            tau_m1 = model.tau_m1(x,u);
            tau_m2 = model.tau_m2(x,u);
            I1 = tau_m1/model.actuator.K1;
            I2 = tau_m2/model.actuator.K2;
            p1_diss = I1^2*model.actuator.R1;
            p2_diss = I2^2*model.actuator.R2;
            %p1_mech = tau_m1*x(5);
            %p2_mech = tau_m2*x(6);
            p1_elec = I1*7.4; %
            p2_elec = I2*7.4; % assume a constant voltage is not right
            power = p1_elec + p2_elec;
            
        end
        
        function [power, p1_elec, p2_elec, p1_diss, p2_diss] = power_elec(model, x, u, varargin)
            parser = inputParser();
            addOptional(parser, 'MotoringOnly', 1);
            parse(parser, varargin{:});
            motoring = parser.Results.MotoringOnly;
            
            tau_m1 = model.tau_m1(x, u);
            tau_m2 = model.tau_m2(x, u);
            I1 = tau_m1/model.actuator.K1;
            I2 = tau_m2/model.actuator.K2;
            p1_diss = I1.^2*model.actuator.R1;
            p2_diss = I2.^2*model.actuator.R2;
            p1_mech = tau_m1.*x(5,:);
            p2_mech = tau_m2.*x(6,:);
            
            
            p1_elec = p1_diss + p1_mech;
            p2_elec = p2_diss + p2_mech;
            
            if motoring
                p1_elec = max(p1_elec,0);
                p2_elec = max(p2_elec,0);
                power = p1_elec + p2_elec;
            else
                power = p1_elec + p2_elec;
            end
            
        end
        


        function power = net_power_load(model,x,u)
            % net output mechanical power on motor level
            % note that motor's power has to be non-negative because energy
            % is not recoverable on motor level
            % Note: don't support vector computation currently
            power_outmech = model.power_load(x,u);
            
            %power_motor1 = sigmf(power_motor1, [10000 0]);
            %power_motor2 = sigmf(power_motor2, [10000 0]);
            
            power = power_outmech - model.power_rege(x,u);
        end
        
        function [power, power_load1, power_load2] = power_load(model,x,u)
            % output mechanical power on motor level
            % note that motor's power has to be non-negative because energy
            % is not recoverable on motor level
            % Note: don't support vector computation currently
%             kappa = model.actuator.Ks;
%             B = model.actuator.B;
%             C = model.actuator.C;
%             r = model.actuator.r;
%             gear = model.actuator.gear;
%             A0 = model.actuator.A0;
%             A = sqrt(B^2+C^2 - 2*B*C*cos(x(3)/gear-x(1)));
%             L = A + r*x(4) - A0; 
            
            tau_s = model.torque_spring(x,u);
            tau_l1 = tau_s/model.actuator.gear;
            tau_l2 = model.actuator.torque_load2(x(1),x(3),x(4));
            
            power_load1 = tau_l1*x(5);
            power_load2 = tau_l2*x(6);
            %if (power_motor1 <= 0)
            %    power_motor1 = 0;
            %end
            %if (power_motor2 <= 0)
            %    power_motor2 = 0;
            %end
            %power_motor1 = sigmf(power_motor1, [10000 0]);
            %power_motor2 = sigmf(power_motor2, [10000 0]);
            power = power_load1 + power_load2;
        end
        
        function [p] = power_link(model, x, u)
            tau_spring = model.torque_spring(x,u);
            p = tau_spring.*x(2,:);
        end
        
        

        %%%% power %%%%
        
        function [ xdot, xdot_x, xdot_u ] = motor_dynamics_2nd(model, x, u)
            % x: theta1 ; theta2; dtheta1; dtheta2;
            % u: u1 u2 - desired positions
            p = model.alpha_servo;
            %tau_l1 = model.torque_spring(x)/model.actuator.gear;
            %tau_l2 = model.actuator.torque_load2(x(1,:),x(3,:),x(4,:));
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
        methods (Static)
            
%             function [xdot, xdot_x, xdot_u] = motor_dynamics_3rd(x,u, p)
%                 z     = zeros(3,1);
%                 Z     = zeros(3,3);
% 
%                 Af1 = [0,     1,     0;
%                     0,     0,     1;
%                     -p(1), -p(2), -p(3)];
% 
%                 Bf1 = [0;
%                     0;
%                     p(1)];
% 
%                 Af2 = [0,     1,     0;
%                     0,     0,     1;
%                     -p(4), -p(5), -p(6)];
%   
%                 Bf2 = [0;
%                    0;
%                    p(4)];
%  
%                 A = [Af1,Z;
%                 Z ,Af2];
%  
%                 B = [Bf1, z;
%                 z ,Bf2];
%     
%                 xdot    = A * x + B * u;  
%                 if nargout > 1
%                 xdot_x  = A;  
%                 xdot_u  = B; 
%                 end
%             end
%             
            
            
            
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