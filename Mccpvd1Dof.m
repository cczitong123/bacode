classdef Mccpvd1Dof < Arm
    %MCCPVD1DOF Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        dimX = 2
        dimU = 3
        M = 0.0036 % inertia
        fv = 0.0077 %viscous friction
        gc = 0 %gravity constant, default 0 for planar
        B = 0.036;
        C = 0.135;
        A0 = 0.099; % C - B
        Ks = 394; %231 % spring constant
        r = 0.015;
        Dm = 0.0339;
        
        mass = 0.18
        com = 0.085
        modelpara
        actuator
        
        param=[]
        %%%% store feedforward command
        %ffu1
        %ffu2
        %%%%
    end
    
    methods
        function obj = Mccpvd1Dof()
            obj = obj@Arm();
            obj.actuator = ActMccpvd();
            obj.modelpara = [];
            obj.modelpara.M = obj.M;
            obj.modelpara.fv = obj.fv;
            obj.modelpara.gc = obj.gc;
            obj.modelpara.B = obj.B;
            obj.modelpara.C = obj.C;
            obj.modelpara.A0 = obj.A0;
            obj.modelpara.Ks = obj.Ks;
            obj.modelpara.r = obj.r;
            obj.modelpara.Dm = obj.Dm;
            obj.modelpara.mass=obj.mass;
            obj.modelpara.com=obj.com;
            
            B = obj.B;
            C = obj.C;
            A0 = obj.C - obj.B; % C - B
            Ks = obj.Ks;
            r = obj.r;
            Dm = obj.Dm;
            
            % q, qdot, th1, th2, th3
            syms F(x1,x2,x3,x4,x5)
            
            F(x1,x2,x3,x4,x5) = Ks*B*C*sin(x3-x1)*(1+ (r*x4 - A0)/sqrt(B^2+C^2-2*B*C*cos(x3-x1)) ) - Dm*x5*x2;
            
            df1 = diff(F,x1);
            df2 = diff(F,x2);
            df3 = diff(F,x3);
            df4 = diff(F,x4);
            df5 = diff(F,x5);
            
            
            df1n = matlabFunction(df1);
            df2n = matlabFunction(df2);
            df3n = matlabFunction(df3);
            df4n = matlabFunction(df4);
            df5n = matlabFunction(df5);
            Fn = matlabFunction(F);
            Ju = [df3,df4,df5];
            Ju = matlabFunction(Ju);
            
            param.Ju = Ju;
            param.Jq = df1n;
            param.Jqdot = df2n;
            
            param.F = Fn;
            param.N = eye(3);
            param.N(3,3)= 1;
            param.inertia= obj.M;
            
            obj.param = param;
        end
        
        function tau = desired_torque(obj, q, qd, qdd)
            %inverse dynamics
            tau = obj.M*qdd + obj.fv*qd - obj.gravity(q);
        end
        
        function f = gravity(obj, q)
            f = -obj.mass*obj.com*obj.gc*sin(q);
            %
        end
        
        function qdd = ffdynamics(obj, q, qd, tau)
            qdd = (tau - obj.gravity(q) - obj.fv*qd)/obj.M;
        end
        
        % SSM dynamics
        % dynamics: u is motor positions
        function xdot = dynamics(obj, x, u)
            qdd = obj.accel(x, u);
            xdot = [x(2);qdd];
        end
        
        % dynamics_v:
        % x: (q, qd, u1, u2, u3)
        % v: motor velocities
        function xdot = dynamics_v(obj, x, v)
            qdd = obj.accel(x(1:2,:), [x(3:4,:); v(3)]);
            xdot = [x(2); qdd; v(1); v(2)];
        end
        
        %         function xdot = dynamics3(obj, x, u3, n)
        %             u1 = obj.ffu1(n);
        %             u2 = obj.ffu2(n);
        %             qdd = obj.accel( x(1:2,:), [u1;u2;u3] );
        %             xdot = [x(2); qdd];
        %         end
        %
        %         function xdot = dynamics2(obj, x, u23, n)
        %             u1 = obj.ffu1(n);
        %             qdd = obj.accel( x(1:2,:), [u1; u23]);
        %             xdot = [ x(2); qdd];
        %         end
        
        function qdd = accel(obj, x, u)
            tau = obj.torque_actuator(x, u);
            %tau = mccp1TauA_mex(obj.modelpara,x,u);
            qdd = (tau - obj.fv*x(2,:) + obj.gravity(x(1,:)))/obj.M;
        end
        
        function tau = torque_actuator(obj, x, u)
            tau = obj.Ks*obj.B*obj.C*sin(u(1,:)-x(1,:)).*(1+ (obj.r*u(2,:) ...
                -obj.A0)./sqrt(obj.B^2+obj.C^2-2*obj.B*obj.C*cos(u(1,:)-x(1,:))) )...
                - obj.Dm*u(3,:).*x(2,:);
        end
        function force = spring_force(obj, x, u)
            
            %sprlength = obj.spring_displacement(q,theta1, theta2);
            force = obj.actuator.spring_force(x(1,:), u(1,:), u(2,:));
        end
        %         function tau = tauA_mex(obj, s, x, u)
        %             tau = mccp1TauA_mex( s, x, u);
        %         end
        % torques
        function tau = torque_spring(obj, x, u)
            tau = obj.actuator.torque_spring(x(1,:), u(1,:), u(2,:));
        end
        function tau = torque_damping(obj, qdot, u3)
            tau = obj.Dm*u3*qdot;
        end
        % motor torques
        function tau_l = tau_l1(obj, x, u)
            tau_l = -obj.torque_spring(x, u);
        end
        
        function tau_l = tau_l2(obj, x, u)
            phi = u(1,:) - x(1,:);
            A = sqrt( obj.C^2 + obj.B^2 - 2*(obj.B*obj.C).*cos( phi ));
            tau_l = -obj.Ks*(A + obj.r*u(2,:) - obj.A0)*obj.r;
            
        end
        
        function tau_m = tau_m1(obj, x, u, yd, ydd)
            % tau_m of servo1 at servo's shaft
            % M ddtheta + D dtheta = tau_m + tau_l
            tau_l = obj.tau_l1(x,u);
            %p = obj.alpha_servo;
            %accel = p^2*(u(1)-x(3)) - 2*p*x(5);
            
            tau_m = - tau_l + obj.actuator.J1*ydd(1,:) + obj.actuator.D1*yd(1,:) ;
            %tau_m = model.actuator.J1*accel + model.actuator.D1*x(5) ;
            %tau_m = tau_l;
        end
        
        function tau_m = tau_m2(obj, x, u, yd, ydd)
            % tau_m of servo2 at servo's shaft
            % yd: motor velocity,
            % ydd: motor's acceleration
            tau_l = obj.tau_l2(x,u);
            %p = obj.alpha_servo;
            %accel = p^2*(u(2)-x(4)) - 2*p*x(6);
            
            tau_m = - tau_l + obj.actuator.D2*yd + obj.actuator.J2*ydd;
            %tau_m = model.actuator.D2*x(6) + model.actuator.J2*accel;
        end
        
        % power
        function p = power_out(obj, x, u)
            % to replace power_link
            tau = obj.torque_spring(x, u);
            p = tau.*x(2,:); % joint velocity
        end
        function out = sqrt_torques(obj, x, u)
            t1 = obj.tau_l1(x,u);
            t2 = obj.tau_l2(x,u);
            out = 5.6611*(t1.^2 + t2.^2); % 5.6611 = R/K_t^2 according to motor specs
            
        end
        % note: before input x, add extra u and udot to x to make up 6 x Nt
        % state matrix
        function p = power_in(obj, x, u)
            %tau_l1 = model.actuator.torque_load1(x(1,:),x(3,:),x(4,:));
            %tau_l2 = model.actuator.torque_load2(x(1,:),x(3,:),x(4,:));
            %p = max(-tau_l1.*x(5,:),0) + max(-tau_l2.*x(6,:),0);
            p1 = obj.power_in1(x, u);
            p2 = obj.power_in2(x, u);
            p1 = max(p1,0);
            p2 = max(p2,0);
            p = p1 + p2;
        end
        function p = power_in1(obj, x, u)
            tau_l1 = obj.actuator.torque_load1(x(1,:),u(1,:),u(2,:));
            p = -tau_l1.*x(5,:); % motor1 velocity
        end
        function p = power_in2(obj, x, u)
            
            tau_l2 = obj.actuator.torque_load2(x(1,:),u(1,:),u(2,:));
            p = -tau_l2.*x(6,:); % motor2 velocity
        end
        function E = energy_spring(obj, x, u)
            l = obj.actuator.spring_displacement(x(1,:),u(1,:),u(2,:));
            E = obj.actuator.Ks*(l.^2)/2;
        end
        function D = evaluate_traj(obj, traj, parameters)
            
        end
        
        function [ traj ] = track_minjerk(obj,q0,qf,tf, varargin)
            parser = inputParser();
            addOptional(parser, 'plot', 0);
            parse(parser, q0, qf, tf, varargin{:});
            
            dt = 0.01;
            simT = 0.1;
            
            [t, q_des, qd_des, qdd_des] = generate_trajectory_jerk(q0, qf, tf-simT, dt);
            
            
            extrat = (tf+dt):dt:simT;
            t = [t,extrat];
            q_des = [q_des, repmat(qf,1,length(extrat))];
            qd_des = [qd_des, zeros(1,length(extrat))];
            qdd_des = [qdd_des, zeros(1,length(extrat))];
            qddd_des = gradient(qdd_des, dt);
            
            
            obj.param.K = [2500 1000 100];
            % initiate states
            x.q = 0;
            x.qd = 0;
            x.qdd = 0;
            x.m1 = 0;
            x.m2 = 0.01;
            x.m3 = 0;
            
            X = zeros(4,length(t));
            
            u = zeros(3, length(t));
            u0 = [x.m1; x.m2; x.m3];
            u(:,1) = u0;
            v = zeros(3, length(t-1));
            vns = zeros(3,1);
            
            BCv = [-7, 7;
                -7, 7;
                -inf, inf];
            BCu = [-pi/2, pi/2;
                0,   pi;
                0,   1];
            
            
            
            for i = 1:length(t)-1 % iterate in time
                % compute desired qdd
                %qdderr(i) = (q_des(i+1)-x.q)/(dt^2);
                %taud(i) = robot.desired_torque(x.q, x.qd, qdderr(i));
                %taudc(i) = robot.desired_torque(x.q, x.qd, x.qdd); % inverse dynamics
                
                %vi = nscontroller_mccpvd(x, taud(i), (taud(i)-taudc(i))/dt, vns, param);
                
                %%%% 1 use desired tau and estimated current tau to determine the change of tau
                %tau_c = robot.torque_actuator([x.q; x.qd], [ x.m1;x.m2;x.m3 ]);
                %taudot_d_i = (tau_d(i)-tau_c)*Kd;
                %%%% 2 use desired position and velocity to determine the change of tau
                
                %%%% 3
                %taudot_d_i = taudot_d(i);
                %%%%
                %vns = [-sqrt(B^2+C^2-2*B*C*cos(x.m1-x.q))*B*C*sin(x.m1-x.q);-r; 0]*10;
                vns = [-(x.m1-qf); -( x.m2-0 ) ; 0];
                %vns = zeros(3,1);
                vi = nscontroller_mccpvd(x, q_des(i),  qd_des(i), qdd_des(i), qddd_des(i), vns, obj.param);
                
                %apply limit on v
                vi = min(max(vi,BCv(:,1)),BCv(:,2));
                %
                uip = u(:,i) + vi*dt;
                % appliy limit on u
                uip = min(max(uip,BCu(:,1)), BCu(:,2));
                % limited v
                vi = (uip - u(:,i))/dt;
                %
                v(:,i) = vi;
                u(:,i+1) = uip;
                % update
                xn = obj.step([x.q; x.qd],u(:,i),dt,dt);
                
                x.q = xn(1);
                x.qd = xn(2);
                
                
                x.m1 = u(1,i+1);
                x.m2 = u(2,i+1);
                x.m3 = u(3,i+1);
                
                X(1,i+1) = x.q;
                X(2,i+1) = x.qd;
                X(3,i+1) = x.m1;
                X(4,i+1) = x.m2;
                
                x.qdd = obj.accel([x.q; x.qd],u(:,i+1));
                
                
                
            end
            
            
            traj.x = X;
            traj.u = u;
            traj.v = v;
            traj.Pin = obj.power_in(traj.x(1:2,:), traj.u);
            traj.Ein = sum(traj.Pin)*dt;
            
            if parser.Results.plot == 1
            figure
            fig = gcf;
            fig.Units='inches';
            fig.Position = [2 1 3.5 6];
            subplot(4,1,4)
            plot(t, u)
            title('Resolved control variables')
            xlabel('Time (s)')
            ylabel('u')
            %legend('u1','u2','u3')
            
            
            subplot(4,1,1)
            title('Position')
            hold on
            plot(t, X(1,:),'-b','LineWidth',1)
            plot(t,q_des,'--r','LineWidth',1)
            legend('result','reference')
            hold off
            ylabel('Position (rad)')
            
            subplot(4,1,2)
            hold on
            plot(t, X(2,:),'-b','LineWidth',1)
            plot(t,qd_des,'--r','LineWidth',1)
            hold off
            title('Velocity')
            ylabel('Velocity (rad/s)')
            
            subplot(4,1,3)
            Xacc = gradient(X(2,:))/0.01;
            hold on
            plot(t, Xacc ,'-b', 'LineWidth',1 )
            plot(t, qdd_des,'--r','LineWidth',1 )
            hold off
            end
        end
        
        function [ traj ] = eidc_track( obj, x0, t, q_des, qd_des, qdd_des, qddd_des, varargin )
            parser = inputParser();
            addOptional(parser, 'plot', 0);
            addOptional(parser, 'm2s', 0.1);
            addOptional(parser, 'qf', NaN);
            addOptional(parser, 'dt', 0.01);
            parse(parser, varargin{:});
            m2s = parser.Results.m2s;
            if isnan(parser.Results.qf)
                qf = q_des(end);
            else
                qf = parser.Results.qf;
            end
            
            dt = parser.Results.dt;
            
            
            obj.param.K = [5000 500 50];
            % initiate states
            x.q = x0(1,1);
            x.qd = x0(2,1);
            x.qdd = x0(3,1);
            x.m1 = x0(4,1);
            x.m2 = x0(5,1);
            x.m3 = x0(6,1);
            
            X = zeros(6,length(t));
            X(:,1)=x0;
            u = zeros(3, length(t));
            u0 = [x.m1; x.m2; x.m3];
            u(:,1) = u0;
            v = zeros(3, length(t)-1);
            %vns = zeros(3,1);
            
            BCv = [-7, 7;
                -7, 7;
                -inf, inf];
            BCu = [-pi/2, pi/2;
                0,   pi;
                0,   1];
            
            
            
            for i = 1:length(t)-1 % iterate in time
                % compute desired qdd
                %qdderr(i) = (q_des(i+1)-x.q)/(dt^2);
                %taud(i) = robot.desired_torque(x.q, x.qd, qdderr(i));
                %taudc(i) = robot.desired_torque(x.q, x.qd, x.qdd); % inverse dynamics
                
                %vi = nscontroller_mccpvd(x, taud(i), (taud(i)-taudc(i))/dt, vns, param);
                
                %%%% 1 use desired tau and estimated current tau to determine the change of tau
                %tau_c = robot.torque_actuator([x.q; x.qd], [ x.m1;x.m2;x.m3 ]);
                %taudot_d_i = (tau_d(i)-tau_c)*Kd;
                %%%% 2 use desired position and velocity to determine the change of tau
                
                %%%% 3
                %taudot_d_i = taudot_d(i);
                %%%%
                %vns = [-sqrt(B^2+C^2-2*B*C*cos(x.m1-x.q))*B*C*sin(x.m1-x.q);-r; 0]*10;
                vns = [-(x.m1-qf); -( x.m2-m2s ) ; 0];
                %vns = zeros(3,1);
                vi = nscontroller_mccpvd(x, q_des(i),  qd_des(i), qdd_des(i), qddd_des(i), vns, obj.param);
                
                %apply limit on v
                vi = min(max(vi,BCv(:,1)),BCv(:,2));
                %
                uip = u(:,i) + vi*dt;
                % appliy limit on u
                uip = min(max(uip,BCu(:,1)), BCu(:,2));
                % limited v
                vi = (uip - u(:,i))/dt;
                %
                v(:,i) = vi;
                u(:,i+1) = uip;
                % update
                xn = obj.step([x.q; x.qd],u(:,i),dt,dt);
                
                x.q = xn(1);
                x.qd = xn(2);
                
                
                x.m1 = u(1,i+1);
                x.m2 = u(2,i+1);
                x.m3 = u(3,i+1);
                
                
                x.qdd = obj.accel([x.q; x.qd],u(:,i+1));
                
                X(1,i+1) = x.q;
                X(2,i+1) = x.qd;
                
                
                X(3,i+1) = x.qdd;
                X(4,i+1) = x.m1;
                X(5,i+1) = x.m2;
                X(6,i+1) = x.m3;
                
                
                
                
            end
            
            X(2,end)=0;
            X(3,end)=0;
            
            traj.X = X;
            traj.x(1:2,:) = X(1:2,:);
            traj.x(3:4,:) = X(4:5,:);
            motorspeed = [v(1:2,:), zeros(2,1)];
           
            traj.x(5:6,:) = motorspeed;
            traj.u = u;
            traj.v = v;
            traj.t = t;
            traj.Pin = obj.power_in(traj.x, traj.u);
            traj.Ein = sum(traj.Pin)*dt;
            traj.Jp =  1000*(qf - traj.x(1,end))^2;
            if parser.Results.plot == 1
            figure
            fig = gcf;
            fig.Units='inches';
            fig.Position = [2 1 3.5 6];
            subplot(4,1,4)
            plot(t, u)
            title('Resolved control variables')
            xlabel('Time (s)')
            ylabel('u')
            %legend('u1','u2','u3')
            
            
            subplot(4,1,1)
            title('Position')
            hold on
            plot(t, X(1,:),'-b','LineWidth',1)
            plot(t,q_des,'--r','LineWidth',1)
            legend('result','reference')
            hold off
            ylabel('Position (rad)')
            
            subplot(4,1,2)
            hold on
            plot(t, X(2,:),'-b','LineWidth',1)
            plot(t,qd_des,'--r','LineWidth',1)
            hold off
            title('Velocity')
            ylabel('Velocity (rad/s)')
            
            subplot(4,1,3)
            Xacc = gradient(X(2,:))/0.01;
            hold on
            plot(t, Xacc ,'-b', 'LineWidth',1 )
            plot(t, qdd_des,'--r','LineWidth',1 )
            hold off
            end
        end
        
    end
    methods (Static)
        
    end
    
    
end

