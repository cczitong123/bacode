classdef Mccpvd2VC < Mccpvd2Dof
    %MCCPVD2VC Maccepa-vd velocity controlled motor
    % x: q1 q2 q1dot q2dot (theta1 theta2)joint1 (theta1 theta2)joint2
    % u: (u1 u2 u3)1 (u1 u2 u3)2
    properties
        name = 'MACCEPA-VD 2 DoF with velocity controlled motors';
        
        alpha_servo = 20; % Fan: fit from data
        
        param
    end
    
    methods
        function obj = Mccpvd2VC()
            %MCCPVD2VC Construct an instance of this class
            %   Detailed explanation goes here
            obj = obj@Mccpvd2Dof();
            obj.actuator1.Ks = 750;
            obj.actuator2.Ks = 750;
            
            B1 = obj.actuator1.B;
            C1 = obj.actuator1.C;
            A01 = obj.actuator1.C - obj.actuator1.B; % C - B
            Ks1 = obj.actuator1.Ks;
            r1 = obj.actuator1.r;
            Dm1 = obj.actuator1.max_damping;
            
            B2 = obj.actuator2.B;
            C2 = obj.actuator2.C;
            A02 = obj.actuator2.C - obj.actuator2.B; % C - B
            Ks2 = obj.actuator2.Ks;
            r2 = obj.actuator2.r;
            Dm2 = obj.actuator2.max_damping;
            
            %      q1, q2, qd1, qd2, m11, m12, u13, m21, m22, u23
            syms F(x1, x2, x3,  x4,   x5,  x6,  x7,  x8,  x9, x10)
            
            F(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10) = [Ks1*B1*C1*sin(x5-x1)*(1+ (r1*x6 - A01)/sqrt(B1^2+C1^2-2*B1*C1*cos(x5-x1)) ) - Dm1*x7*x3;
                Ks2*B2*C2*sin(x8-x2)*(1+ (r2*x9 - A02)/sqrt(B2^2+C2^2-2*B2*C2*cos(x8-x2)) ) - Dm2*x10*x4 ];
            
            df1 = diff(F,x1);
            df2 = diff(F,x2);
            df3 = diff(F,x3);
            df4 = diff(F,x4);
            df5 = diff(F,x5);
            df6 = diff(F,x6);
            df7 = diff(F,x7);
            df8 = diff(F,x8);
            df9 = diff(F,x9);
            df10 = diff(F,x10);
            
            
            Fn = matlabFunction(F);
            
            Ju = [df5,df6,df7,df8,df9,df10];
            Ju = matlabFunction(Ju);
            Jq = [df1,df2];
            Jq = matlabFunction(Jq);
            Jqdot = [df3, df4];
            Jqdot = matlabFunction(Jqdot);
            
            obj.param.Ju = Ju;
            obj.param.Jq = Jq;
            obj.param.Jqdot = Jqdot;
            
            obj.param.F = Fn;
            obj.param.N = diag([1 1 0.1 1 1 0.1]);
            
        end
        
        function [ traj ] = track_joint_space(obj, x0, q_des, varargin)
            
            parser = inputParser();
            
            %addRequired(parser, 'robot');
            addRequired(parser, 'x0');
            addRequired(parser, 'q_des');
            addOptional(parser, 'qd_des', NaN);
            addOptional(parser, 'qdd_des', NaN);
            addOptional(parser, 'qddd_des', NaN);
            addOptional(parser, 'dt', 0.01);
            addOptional(parser, 't', NaN);
            addOptional(parser, 'qf', NaN);
            addOptional(parser, 'Textra', 0);
            addOptional(parser, 'K', [2500 1000 100]);
            addOptional(parser, 'plot', 0);
            addOptional(parser, 'm2s',[0;0]);
            addOptional(parser, 'wns', 1); % null-space control weight
            parse(parser, x0, q_des, varargin{:});
            
            
            obj.param.K = parser.Results.K;
            dt = parser.Results.dt;
            q_des = parser.Results.q_des;
            m2s = parser.Results.m2s;
            wns = parser.Results.wns;
            if isnan(parser.Results.qd_des)
                qd_des = gradient(q_des)/dt;
            else
                qd_des = parser.Results.qd_des;
            end
            if isnan(parser.Results.qdd_des)
                qdd_des = gradient(qd_des)/dt;
            else
                qdd_des = parser.Results.qdd_des;
            end
            if isnan(parser.Results.qddd_des)
                qddd_des = gradient(qdd_des)/dt;
            else
                qddd_des = parser.Results.qddd_des;
            end
            
            if isnan(parser.Results.qf)
                qf = q_des(1:2,end);
            else
                qf = parser.Results.qf;
            end
            if isnan(parser.Results.t)
                t = 0:dt:( (size(q_des,2)-1 )*dt+ parser.Results.Textra);
            else
                t = parser.Results.t;
            end
            
            % initiate states
            %x.q = [0;0];
            %x.qd = 0;
            %x.qdd = 0;
            %x.m1 = 0;
            %x.m2 = 0.01;
            %x.m3 = 0;
            x = zeros(8,length(t));
            x(:,1) = x0;
            
            v = zeros(6, length(t)-1);
            
            u = zeros(6, length(t));
            u(1,1) = x0(5);
            u(2,1) = x0(6);
            u(4:5,1) = x0(7:8);
            
            xx.q = x(1:2,1);
            xx.qd = x(3:4,1);
            uu = v(:,1); uu(3) = u(3,1); uu(6) = u(6,1);
            xx.qdd = obj.qddot(x(:,1), uu);
            traj.qdd(:,1)=xx.qdd;
            
            
            %u0 = [x.m1; x.m2; x.m3];
            
            %u(:,1) = u0;
            
            %vns = zeros(6,1);
            
            BCv = [-7, 7;
                -7, 7;
                -inf, inf;
                -7, 7;
                -7, 7;
                -inf, inf];
            
            BCu = [-pi/2, pi/2;
                0,   pi;
                0,   1;
                -pi/2, pi/2;
                0,   pi;
                0,   1];
            
            
            xs = zeros(10,1);
            
            
            
            for i = 1:length(t)-1 % iterate in time
                xs(1:4) = x( 1:4 , i );
                xs(5:10) = u(:, i) ;
                %%%%
                %vns = [-sqrt(B^2+C^2-2*B*C*cos(x.m1-x.q))*B*C*sin(x.m1-x.q);-r; 0]*10;
                vns = wns*[ -( x(5,i) - qf(1) ); -( x(6,i) - m2s(1) ) ; 0; -( x(7,i) - qf(2) ); -( x(8,i) - m2s(2) ); 0 ];
                %vns = zeros(3,1);
                vi = obj.nscontroller( xs, xx, q_des(:,i),  qd_des(:,i), qdd_des(:,i), qddd_des(:,i), vns, obj.param);
                
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
                
                % update
                uu = vi; uu(3) = u(3,i); uu(6) = u(6,i);
                xn = obj.step(x(:,i) , uu, dt, dt);
                
                u(:,i+1) = uip;
                x(:,i+1) = xn;
                
                xx.q = xn(1:2);
                xx.qd = xn(3:4);
                uu = vi; uu(3) = u(3,i+1); uu(6) = u(6,i+1);
                xx.qdd = obj.qddot( xn, uu );
                traj.qdd(:,i+1) = xx.qdd;
                %x.qdd = robot.accel([x.q; x.qd],u(:,i+1));
                
                
                
            end
            traj.x = x;
            traj.u = v;
            traj.u(3,:) = u(3,1:end-1) ;
            traj.u(6,:) = u(6,1:end-1) ;
            traj.Pin = obj.power_in(traj.x(:,1:end-1),traj.u);
            traj.Ein = sum(traj.Pin)*dt;
            traj.Pelec = obj.power_elec(traj.x(:,1:end-1),traj.u);
            traj.Eelec = sum(traj.Pelec)*dt;
            traj.q = traj.x(1:2,:);
            traj.qd = traj.x(3:4,:);
            traj.t = t;
            traj.q_des= q_des;
            traj.qd_des=qd_des;
            traj.qdd_des=qdd_des;
            %traj.qdd = gradient(  );
            %traj.v = v;
        end
        
        function [ traj ] = reach_minjerk_os(obj, x0, y0, yf, tf, varargin)
                
            parser = inputParser();
            
            %addRequired(parser, 'robot');
            addRequired(parser, 'x0');
            addRequired(parser, 'y0');
            addRequired(parser, 'yf');
            addRequired(parser, 'tf');
            %addOptional(parser, 'yddd_des', NaN);
            addOptional(parser, 'dt', 0.01);
            addOptional(parser, 't', NaN);
            addOptional(parser, 'Textra', 0.1);
            addOptional(parser, 'K', [2500 1000 100]);
            addOptional(parser, 'plot', 0);
            addOptional(parser, 'm2s',[0.1;0.1]);
            addOptional(parser, 'animate' , 0);
            addOptional(parser, 'save_animate' , 0);
            parse(parser, x0, y0, yf, tf, varargin{:});
            K  = parser.Results.K;
            x0 = parser.Results.x0;
            tf = parser.Results.tf;
            dt = parser.Results.dt;
            y0 = parser.Results.y0;
            yf = parser.Results.yf;
            
            [ t, y_des, yd_des, ydd_des ] = generate_trajectory_jerk(y0, yf, tf, dt);
            
            simT = tf + parser.Results.Textra;
            extrat = (tf+dt):dt:simT;
            t = [t,extrat];
            y_des = [y_des, repmat(yf,1,length(extrat))];
            yd_des = [yd_des, zeros(2,length(extrat))];
            ydd_des = [ydd_des, zeros(2,length(extrat))];
            yddd_des = gradient(ydd_des, dt);
            traj = obj.track_os_space( x0, y_des, 'yd_des', yd_des , 'ydd_des', ydd_des , 'yddd_des', yddd_des,'yf',yf,'t',t,'m2s',parser.Results.m2s,'K',K );
            
            if parser.Results.plot == 1
                figure
                fig = gcf;
                set(fig,'color','white')
                fig.Units='inches';
                fig.Position=[5 2 3.5 6];
                subplot(3,1,1)
                plot(traj.t, traj.x(1,:))
                hold on
                plot(traj.t, traj.x(2,:))
                title('Joints')
                ylabel('q_{1,2} [rad]')
                legend('Joint 1','Joint 2')
                subplot(3,1,2)
                plot(traj.t, traj.x(5,:))
                hold on
                plot(traj.t, traj.x(7,:))
                title('EP motors')
                ylabel('[rad]')
                
                subplot(3,1,3)
                plot(traj.t, traj.x(6,:) )
                hold on
                plot(traj.t, traj.x(8,:) )
                xlabel('t [s]')
                ylabel('[rad]')
                title('Stiffness motors')
            end
            
            if parser.Results.animate == 1
                frames = obj.animate(traj.x, 0.01, 'goal', yf);
                if parser.Results.save_animate == 1
                Arm2Dof.savegif(frames, 'Mccpvd2VC_eidc_minjerkreach.gif')
                end
            end
        end
        
        function [ traj ] = track_os_space(obj, x0, y_des, varargin)
            
            parser = inputParser();
            
            %addRequired(parser, 'robot');
            addRequired(parser, 'x0');
            addRequired(parser, 'y_des');
            addOptional(parser, 'yd_des', NaN);
            addOptional(parser, 'ydd_des', NaN);
            addOptional(parser, 'yddd_des', NaN);
            addOptional(parser, 'dt', 0.01);
            addOptional(parser, 't', NaN);
            addOptional(parser, 'yf', NaN);
            
            addOptional(parser, 'K', [100 0 0]);
            addOptional(parser, 'plot', 0);
            addOptional(parser, 'm2s', [0;0])
            parse(parser, x0, y_des, varargin{:});
            
            m2s = parser.Results.m2s;
            obj.param.K = parser.Results.K;
            dt = parser.Results.dt;
            y_des = parser.Results.y_des;
            if isnan(parser.Results.yd_des)
                yd_des = gradient(y_des)/dt;
            else
                yd_des = parser.Results.yd_des;
            end
            if isnan(parser.Results.ydd_des)
                ydd_des = gradient(yd_des)/dt;
            else
                ydd_des = parser.Results.ydd_des;
            end
            if isnan(parser.Results.yddd_des)
                yddd_des = gradient(ydd_des)/dt;
            else
                yddd_des = parser.Results.yddd_des;
            end
            
            if isnan(parser.Results.yf)
                yf = y_des(1:2,end);
            else
                yf = parser.Results.yf;
            end
            if isnan(parser.Results.t)
                t = 0:dt:( (size(y_des,2)-1 )*dt+ parser.Results.Textra);
            else
                t = parser.Results.t;
            end
            
            % initiate states
            %x.q = [0;0];
            %x.qd = 0;
            %x.qdd = 0;
            %x.m1 = 0;
            %x.m2 = 0.01;
            %x.m3 = 0;
            x = zeros(8,length(t));
            x(:,1) = x0;
            
            v = zeros(6, length(t)-1);
            
            u = zeros(6, length(t));
            u(1,1) = x0(5);
            u(2,1) = x0(6);
            u(4:5,1) = x0(7:8);
            
            xx.q = x(1:2,1);
            xx.qd = x(3:4,1);
            uu = v(:,1); uu(3) = u(3,1); uu(6) = u(6,1);
            xx.qdd = obj.qddot(x(:,1), uu);
            xx.qddd = [0;0];
            xx.y = obj.endpoint(xx.q);
            xx.yd = Arm2Dof.jacobian(xx.q,obj.L)*xx.qd;
            
            
            %u0 = [x.m1; x.m2; x.m3];
            
            %u(:,1) = u0;
            
            %vns = zeros(6,1);
            
            BCv = [-7, 7;
                -7, 7;
                -inf, inf;
                -7, 7;
                -7, 7;
                -inf, inf];
            
            BCu = [-pi/2, pi/2;
                0,   pi;
                0,   1;
                -pi/2, pi/2;
                0,   pi;
                0,   1];
            
            
            xs = zeros(10,1);
            
            
            
            for i = 1:length(t)-1 % iterate in time
                xs(1:4) = x( 1:4 , i );
                xs(5:10) = u(:, i) ;
                %%%%
                %vns = [-sqrt(B^2+C^2-2*B*C*cos(x.m1-x.q))*B*C*sin(x.m1-x.q);-r; 0]*10;
                vns = [ 0; -( x(6,i) - m2s(1) ) ; 0; 0; -( x(8,i) - m2s(2) ); 0 ];
                %vns = zeros(3,1);
                vi = obj.eidcontroller_os( xs, xx, y_des(:,i),  yd_des(:,i), ydd_des(:,i),yddd_des(:,i) , vns, obj.param);
                
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
                
                % update
                uu = vi; uu(3) = u(3,i); uu(6) = u(6,i);
                xn = obj.step(x(:,i) , uu, dt, dt);
                
                u(:,i+1) = uip;
                x(:,i+1) = xn;
                
                xx.q = xn(1:2);
                xx.qd = xn(3:4);
                uu = vi; uu(3) = u(3,i+1); uu(6) = u(6,i+1);
                oldqdd = xx.qdd;
                xx.qdd = obj.qddot( xn, uu );
                xx.qddd = (xx.qdd-oldqdd)/dt;
                %x.qdd = robot.accel([x.q; x.qd],u(:,i+1));
                xx.y = obj.endpoint(xx.q);
                xx.yd = Arm2Dof.jacobian(xx.q,obj.L)*xx.qd;
                
                
            end
            traj.x = x;
            traj.u = v;
            traj.u(3,:) = u(3,1:end-1) ;
            traj.u(6,:) = u(6,1:end-1) ;
            traj.Pin = obj.power_in(traj.x(:,1:end-1),traj.u);
            traj.Ein = sum(traj.Pin)*dt;
            traj.Pelec = obj.power_elec(traj.x(:,1:end-1),traj.u);
            traj.Eelec = sum(traj.Pelec)*dt;
            traj.t = t;
            traj.y_des = y_des;
            traj.yd_des = yd_des;
            traj.ydd_des = ydd_des;
            traj.yddd_des = yddd_des;
            traj.y = obj.endpoint(traj.x);
            traj.yd(1,:) = gradient( traj.y(1,:),dt );
            traj.yd(2,:) = gradient( traj.y(2,:),dt );
            traj.ydd(1,:) = gradient( traj.yd(1,:),dt );
            traj.ydd(2,:) = gradient( traj.yd(2,:),dt );
            
            %traj.v = v;
        end
        function [ v ] = nscontroller(obj, x, xx, q_des,  qd_des, qdd_des, qddd_des, v1, param)
            q = xx.q;
            qdot = xx.qd;
            qddot = xx.qdd;
            
            dm = obj.M(2)*obj.L(1)*obj.Lg(2)*sin(q(2))*qdot(2);
            dM = [-2*dm, -dm;
                -dm,   0 ];
            CMatrix = Arm2Dof.compute_C(q, qdot,obj);
            dC = obj.M(2)*obj.L(1)*obj.Lg(2)*sin(q(2))*[ -2*qddot(2), -qddot(2); qddot(1), 0  ] + ...
                obj.M(2)*obj.L(1)*obj.Lg(2)*cos(q(2))*[ -2*qdot(2), -qdot(2); qdot(1), 0 ];
            dG = obj.g*[ obj.M(2)*obj.L(1)*sin(q(1)+q(2))*(qdot(1)+qdot(2)) + obj.M(2)*obj.Lg(2)*sin(q(1)+q(2))*(qdot(1)+qdot(2)) ; ...
                obj.M(2)*obj.Lg(2)*sin(q(1)+q(2))*(qdot(1)+qdot(2))  ];
            MM = Arm2Dof.compute_Mq(q, obj);
            %b = taudot_d - param.Jq(q,qdot,u1,u2,u3)*qdot - param.Jqdot(q,qdot,u1,u2,u3)*qddot + param.Kp*(tau_d - param.F(q, qdot, u1, u2, u3));
            b = MM*( qddd_des + param.K(3)*(qdd_des - qddot)+param.K(2)*(qd_des - qdot)+param.K(1)*(q_des-q) )...
                - param.Jq(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10))*qdot - param.Jqdot(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9) ,x(10))*qddot ...
                + dM*qddot + CMatrix*qddot + dC*qdot + dG;
            
            J = param.Ju(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10));
            v = nscontroller(J,b,v1,param.N);
        end
        function [ v ] = eidcontroller_os(obj, x, xx, y_des,  yd_des, ydd_des, yddd_des, v1, param)
            q = xx.q;
            qd = xx.qd;
            qdd = xx.qdd;
            %qddd = xx.qddd;
            y = xx.y;
            yd = xx.yd;
            
            Jk = Arm2Dof.jacobian(q, obj.L);
            dJ = Arm2Dof.compute_dJ(q, qd, obj.L);
            ddJ = Arm2Dof.compute_ddJ(q, qdd, obj.L);
            dm = obj.M(2)*obj.L(1)*obj.Lg(2)*sin(q(2))*qd(2);
            dM = [-2*dm, -dm;
                -dm,   0 ]; % dM/dt
            CM = Arm2Dof.compute_C(q, qd,obj); % Col Matrix
            dC = obj.M(2)*obj.L(1)*obj.Lg(2)*sin(q(2))*[ -2*qdd(2), -qdd(2); qdd(1), 0  ] + ...
                obj.M(2)*obj.L(1)*obj.Lg(2)*cos(q(2))*[ -2*qd(2), -qd(2); qd(1), 0 ];
            dG = obj.g*[ obj.M(2)*obj.L(1)*sin(q(1)+q(2))*(qd(1)+qd(2)) + obj.M(2)*obj.Lg(2)*sin(q(1)+q(2))*(qd(1)+qd(2)) ; ...
                obj.M(2)*obj.Lg(2)*sin(q(1)+q(2))*(qd(1)+qd(2))  ];
            MM = Arm2Dof.compute_Mq(q, obj); % Inertia Matrix
            
            Jqc  = param.Jq(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10));
            Jqdot =  param.Jqdot(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9) ,x(10))  ;
            AA = Jk/MM;
            
            b =  yddd_des + param.K(3)*(ydd_des - ydd_des) + param.K(2)*(yd_des - yd) + param.K(1)*(y_des - y) ...
                - 2*dJ*qdd - ddJ*qd - AA*( Jqc*qd + Jqdot*qdd - dM*qdd - CM*qdd - dC*qd - dG) ;
            
            J = param.Ju(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10));
            A = AA*J ;
            v = nscontroller(A,b,v1,param.N);
        end
        function xdot = dynamics(obj, x, u)
            % uu: (m1 m2 u3)_1 , (m1 m2 u3)_2
            %uu = [x(5:6,:); u(3,:); x(9:10,:); u(6,:)];
            % joint accel
            qddot = obj.qddot( x, u);
            % first 4 of xdot
            xdot1 = [x(3:4,:); qddot];
            % last 8 of xdot
            xdot2 = u(1:2,:);
            
            xdot3 = u(3:4,:);
            
            xdot = [xdot1; xdot2; xdot3];
            
        end
        
        function [xdot, xdot_x, xdot_u] = dynamics_with_jacobian(obj, x, u)
            % uu: (m1 m2 u3)_1 , (m1 m2 u3)_2
            %uu = [x(5:6,:); u(3,:); x(9:10,:); u(6,:)];
            % joint accel
            %qddot = model.qddot( x, u);
            % first 4 of xdot
            %xdot1 = [x(3:4); qddot];
            % last 8 of xdot
            %[xdot2, xdot2_x, xdot2_u] = model.motor_dynamics_2nd(x(5:8,:), u(1:2,:));
            
            %[xdot3, xdot3_x, xdot3_u] = model.motor_dynamics_2nd(x(9:12,:), u(4:5,:));
            
            %xdot = [xdot1; xdot2; xdot3];
            xdot = obj.dynamics(x,u);
            if nargout > 1
                fx = @(x)obj.dynamics(x, u);
                dfdx = get_jacobian_fd(fx, x);
                if iscolumn(dfdx), dfdx = dfdx'; end
                % Compute derivative of acceleration w.r.t. u
                %daccdu = zeros(1,model.dimU);
                fu = @(u)obj.dynamics(x, u);
                dfdu = get_jacobian_fd(fu, u);
                if iscolumn(dfdu), dfdu = dfdu'; end
                xdot_x = dfdx;
                
                xdot_u = dfdu;
            end
        end
        function acc = qddot(obj, x, u)
            uu = [x(5:6,:); u(3,:); x(7:8,:); u(6,:)];
            q = x(1:2,:);
            qdot = x(3:4,:);
            
            tau = obj.tau(q, qdot, uu);
            acc = Arm2Dof.compute_qddot(q, qdot, tau, obj);
        end
        function x = endpoint(obj, q)
            x = Arm2Dof.endpoint(q, obj.L);
        end
        
        function J = jacobian(obj, q)
            J = Arm2Dof.jacobian(q, obj.L);
        end
        
        function k = stiffness(obj, x)
            k1 = obj.actuator1.stiffness(x(1,:),x(5,:),x(6,:));
            k2 = obj.actuator2.stiffness(x(2,:),x(7,:),x(8,:));
            k = [k1; k2];
        end
        
        %variable damping
        function d = damping(obj, u)
            % duty_circle: 0-1
            % linear on duty-circle
            d1 = obj.actuator1.damping(u(3,:));
            d2 = obj.actuator2.damping(u(6,:));
            d = [d1;d2];
        end
        
        function dr = damping_ratio(obj, x, u)
            b = [obj.viscous_friction;obj.viscous_friction] + obj.damping(u);
            k = obj.stiffness(x);
            dr = b./(2*sqrt(k*obj.inertia));
        end
        
        function [ p, P ] = power_in(obj, x, u)
            tau_11 = obj.actuator1.torque_load1(x(1,:),x(5,:),x(6,:));
            tau_12 = obj.actuator1.torque_load2(x(1,:),x(5,:),x(6,:));
            tau_21 = obj.actuator2.torque_load1(x(2,:),x(7,:),x(8,:));
            tau_22 = obj.actuator2.torque_load2(x(2,:),x(7,:),x(8,:));
            p1 = max(-tau_11.*u(1,:), 0);
            p2 = max(-tau_12.*u(2,:), 0);
            p3 = max(-tau_21.*u(4,:), 0);
            p4 = max(-tau_22.*u(5,:), 0);
            P = [p1;p2;p3;p4];
            p = p1 + p2 + p3 + p4;
        end
        
        function [power, P] = power_elec(obj, x, u, varargin)
            parser = inputParser();
            addOptional(parser, 'MotoringOnly', 1);
            parse(parser, varargin{:});
            motoring = parser.Results.MotoringOnly;
            
            tau_11 = obj.actuator1.torque_load1(x(1,:),x(5,:),x(6,:));
            tau_12 = obj.actuator1.torque_load2(x(1,:),x(5,:),x(6,:));
            tau_21 = obj.actuator2.torque_load1(x(2,:),x(7,:),x(8,:));
            tau_22 = obj.actuator2.torque_load2(x(2,:),x(7,:),x(8,:));
            %tau_m1 = obj.tau_m1(x, u);
            %tau_m2 = obj.tau_m2(x, u);
            tau_m11 = - tau_11 + obj.actuator1.D1*u(1,:) ;
            tau_m12 = - tau_12 + obj.actuator1.D2*u(2,:) ;
            tau_m21 = - tau_21 + obj.actuator2.D1*u(4,:) ;
            tau_m22 = - tau_22 + obj.actuator2.D2*u(5,:) ;
            
            
            I1 = tau_m11/obj.actuator1.K1;
            I2 = tau_m12/obj.actuator1.K2;
            I3 = tau_m21/obj.actuator2.K1;
            I4 = tau_m22/obj.actuator2.K2;
            
            p1_diss = I1.^2*obj.actuator1.R1;
            p2_diss = I2.^2*obj.actuator1.R2;
            p3_diss = I3.^2*obj.actuator2.R1;
            p4_diss = I4.^2*obj.actuator2.R2;
            
            p1_mech = max(tau_m11.*u(1,:),0);
            p2_mech = max(tau_m12.*u(2,:),0);
            
            p3_mech = max(tau_m21.*u(4,:),0);
            p4_mech = max(tau_m22.*u(5,:),0);
            
            p1_elec = p1_diss + p1_mech;
            p2_elec = p2_diss + p2_mech;
            p3_elec = p3_diss + p3_mech;
            p4_elec = p4_diss + p4_mech;
            %if motoring
            %    p1_elec = max(p1_elec,0);
            %    p2_elec = max(p2_elec,0);
            %    power = p1_elec + p2_elec;
            %else
            %    power = p1_elec + p2_elec;
            %end
            power = p1_elec + p2_elec + p3_elec + p4_elec ;
            P = [p1_elec ; p2_elec ; p3_elec ; p4_elec];
        end
        
        function [l, l_x, l_xx, l_u, l_uu, l_ux] = costfnt_reach(obj, x,u,t, costpara)
            %Rtf  = diag([w_tf, 0, 0, 0, 0, 0]);
            Hf = costpara.Hf;
            Hr = costpara.Hr;
            xf = costpara.xf;
            if (isnan(u))
                % final cost
                %fh = @(x)obj.costf_reach(x);
                
                l = ((x-xf)')*Hf*(x-xf)/2 ; %+ Hf(5)*u(1)^2/2 + Hf(6)*u(2)^2/2 + Hf(7)*u(4)^2/2 + Hf(8)*u(5)^2/2;
                if nargout > 1
                    %l_x = get_jacobian_fd(fh, x);
                    %l_xx = get_hessian_fd(fh, x);
                    l_x = ((x-xf)')*Hf ;
                    l_xx = Hf;
                end
            else
                % running cost
                %para = [];
                %para.w_t = w_t;
                %para.w_e = w_e;
                %fl = @(x,u,t) obj.costr_reach(x,u,t);
                %l = fl(x,u,t);
                l = u'*Hr*u;
                if nargout > 1
                    % finite difference
                    %flJ=@(x,u,t)J_cost_fd ( fl, x, u, t );
                    %[l_x ,l_u      ] = flJ ( x, u, t );
                    %flH =@(x,u,t)H_cost_fd  ( flJ, x, u, t );
                    %[l_xx,l_uu,l_ux] = flH  ( x, u, t );
                    
                    l_x = [0 0 0 0 0 0 0 0];
                    l_u = u'*Hr;
                    l_xx = zeros(8,8);
                    l_ux = zeros(6,8);
                    l_uu = Hr;
                    
                end
            end
        end
        
        function [l, l_x, l_xx, l_u, l_uu, l_ux] = costfnt_reachOs(obj, x,u,t, costpara)
            %Rtf  = diag([w_tf, 0, 0, 0, 0, 0]);
            %Hy = costpara.Hy;
            Hr = costpara.Hr;
            %yf = costpara.yf;
            if (isnan(u))
                % final cost
                fh = @(x)obj.costf_reachOs(x, costpara);
                l = fh(x);
                %y = obj.endpoint( x(1:2,:) );
                %l = ((y-yf)')*Hy*(y-yf)/2 ; %+ Hf(5)*u(1)^2/2 + Hf(6)*u(2)^2/2 + Hf(7)*u(4)^2/2 + Hf(8)*u(5)^2/2;
                if nargout > 1
                    l_x = get_jacobian_fd(fh, x);
                    l_xx = get_hessian_fd(fh, x);
                    %l_x = [((y-yf)')*Hy*(obj.jacobian(x(1:2,:)))',0, 0, 0, 0,] ;
                    %l_xx = Hf;
                end
            else
                % running cost
                %para = [];
                %para.w_t = w_t;
                %para.w_e = w_e;
                %fl = @(x,u,t) obj.costr_reach(x,u,t);
                %l = fl(x,u,t);
                l = u'*Hr*u;
                if nargout > 1
                    % finite difference
                    %flJ=@(x,u,t)J_cost_fd ( fl, x, u, t );
                    %[l_x ,l_u      ] = flJ ( x, u, t );
                    %flH =@(x,u,t)H_cost_fd  ( flJ, x, u, t );
                    %[l_xx,l_uu,l_ux] = flH  ( x, u, t );
                    
                    l_x = [0 0 0 0 0 0 0 0];
                    l_u = u'*Hr;
                    l_xx = zeros(8,8);
                    l_ux = zeros(6,8);
                    l_uu = Hr;
                    
                end
            end
        end
        
        function [ traj ] = reach_ilqr(obj, x0, qf, varargin)
            parser = inputParser();
            
            %addRequired(parser, 'robot');
            addRequired(parser, 'x0');
            addRequired(parser, 'qf');
            addOptional(parser, 'tf', 0.8);
            addOptional(parser, 'dt', 0.01);
            addOptional(parser, 'm2f', [0,0]);
            addOptional(parser, 'sfactor', 1);
            %addOptional(parser, 'ConstraintType', 1);
            addOptional(parser, 'plot', 0);
            addOptional(parser, 'animate', 0);
            addOptional(parser, 'saveanimate', 0);
            addOptional(parser, 'os', 0); % reach in opertional space
            parse(parser, x0, qf, varargin{:});
            
            dt = parser.Results.dt;
            tf = parser.Results.tf;
            m2f = parser.Results.m2f;
            
            t = 0:dt:tf;
            Nt = tf/dt + 1;
            
            
            f = @(x,u)obj.dynamics_with_jacobian( x, u);
            costpara = [];
            
            if parser.Results.os == 0
                
                costpara.Hf = 1e3*diag([1 1 1 1 1 1 1 1]);
                %costpara.Hf = 1000*diag([1 1 1 1 0 0 0 0]);
                costpara.Hr = 1*diag([1 1 0.001 1 1 0.001])*parser.Results.sfactor;
                costpara.qf = qf;
                costpara.xf = [ qf(1);qf(2);0;0; qf(1);m2f(1);qf(2);m2f(2) ];
                costpara.tf = tf;
                
                %costpara.L = robot.L ; % link lengths
                %%%%
                j = @(x,u,t)obj.costfnt_reach(x, u, t, costpara);
                
            else
                costpara.yf = qf;
                costpara.Hy = 1e0*diag([1 1]);
                costpara.Hr = 1e-4*diag([1 1 0.001 1 1 0.001])*parser.Results.sfactor;
                costpara.tf = tf;
                j = @(x,u,t)obj.costfnt_reachOs(x, u, t, costpara);
            end
            
            %%%% Define parameters of optimisation
            opt_param = [];
            
            opt_param.lambda_init = 1;
            opt_param.lambda_max  = 1e10;
            opt_param.iter_max = 100;
            opt_param.online_plotting = 0;
            opt_param.online_printing = 2;
            opt_param.dcost_converge = 1e-6;
            opt_param.solver = 'rk4';
            
            opt_param.T = tf;
            %%%
            opt_param.umax = [ 7 ; 7 ; 1 ;  7;  7; 1] ;
            opt_param.umin = [-7; -7 ; 0 ; -7; -7; 0] ;
            % u0 can be full command sequence or just initial point
            %u0 = [qf(1); m2min(1)+0.1; 0; qf(2) ;m2min(2)+0.1 ; 0];
            if parser.Results.os==0
                u0 = [ qf(1)-x0(1); 0 ; 0; qf(2)-x0(2); 0; 0];
            else
                u0= [ randn(1); randn(1) ; 0; randn(1); randn(1); 0];
                
            end
            %%%%
            tic
            %traj = ILQRController.ilqr(f, j, dt, Nt, x0, u0, opt_param);
            
            traj = iLQRv1_mccp2vc(f, j, dt, Nt, x0, u0, opt_param);
            traj.t = t;
            %traj = val_traj_mccpvd(robot, task, traj);
            %if nargin == 0
            %    plot_traj_mccpvd1(traj);
            %end
            
            toc
            disp('ILQR done')
            
            traj.stiffness = obj.stiffness( traj.x );
            traj.damping = obj.damping(traj.u);
            sdt = 0.001;
            tsim = 0:sdt:traj.t(end); % simulation timestamp
            usim = scale_controlSeq(traj.u, traj.t, tsim);
            psim.solver = 'rk4';
            psim.dt = sdt;
            %f = @(x,u) robot.dynamics(x,u);
            xsim = simulate_feedforward(f, traj.x(:,1), usim, psim);
            traj.tsim = tsim;
            traj.usim = usim;
            traj.xsim = xsim;
            traj.qf = qf;
            %if isnan(traj.qf)
            %    sinfo = stepinfo(traj.xsim(1:2,:), [traj.tsim;traj.tsim], 'SettlingTimeShreshold', 0.02);
            %else
            %    sinfo = stepinfo(traj.xsim(1:2,:), [traj.tsim;traj.tsim], traj.qf, 'SettlingTimeShreshold', 0.02);
            %end
            %traj.SettlingTime = sinfo.SettlingTime;
            if parser.Results.os == 0
                traj.Jp = (traj.x(:,end)-costpara.xf)'*costpara.Hf*(traj.x(:,end) - costpara.xf)/2;
            else
                traj.y = obj.endpoint(traj.x);
                traj.Jp = (traj.y(:,end)-costpara.yf)'*costpara.Hy*(traj.y(:,end)-costpara.yf)/2;
            end
            
            [traj.Pin, traj.PinM] = obj.power_in( traj.xsim(:,1:end-1) , traj.usim);
            [traj.Pelec, traj.PelecM] = obj.power_elec( traj.xsim(:,1:end-1) , traj.usim);
            traj.Ein = sum(traj.Pin)*sdt;
            traj.Eelec = sum(traj.Pelec)*sdt;
            
            
            if parser.Results.plot == 1
                figure
                fig = gcf;
                set(fig,'color','white')
                fig.Units='inches';
                fig.Position=[5 2 3.5 6];
                subplot(3,1,1)
                plot(traj.t, traj.x(1,:))
                hold on
                plot(traj.t, traj.x(2,:))
                title('Joints')
                ylabel('q_{1,2} [rad]')
                legend('Joint 1','Joint 2')
                subplot(3,1,2)
                plot(traj.t, traj.x(5,:))
                hold on
                plot(traj.t, traj.x(7,:))
                title('EP motors')
                ylabel('[rad]')
                
                subplot(3,1,3)
                plot(traj.t, traj.x(6,:) )
                hold on
                plot(traj.t, traj.x(8,:) )
                xlabel('t [s]')
                ylabel('[rad]')
                title('Stiffness motors')
            end
            
            if parser.Results.animate == 1
                frames = obj.animate(traj.x, 0.01, 'goal', obj.endpoint(qf));
                Arm2Dof.savegif(frames, 'Mccpvd2VC_reaching.gif')
            end
        end
        
        
        function [ traj ] = evaluate_traj(obj, traj )
            
        end
        function [l] = costf_reachOs(obj,x,costpara)
            Hy = costpara.Hy;
            yf = costpara.yf;
            y = obj.endpoint( x(1:2,:) );
            l =  ((y-yf)')*Hy*(y-yf)/2 ;
        end
        function [ duration]= determine_duration(obj, y0, yf)
            
        end
    end
end

