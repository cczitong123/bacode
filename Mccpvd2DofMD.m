classdef Mccpvd2DofMD < Mccpvd2Dof
    % Mccpvd2 with 2nd order dynamics
    % x: q1 q2 qd1 qd2 (theta1 theta2 dtheta1 dtheta2)_1 (theta1 theta2
    % dtheta1 dtheta2)_2
    % dimX = 12
    % dimU = 6
    
    properties
        name = 'MACCEPA-VD 2 DoF with 2nd order motor dynamics';
        
        alpha_servo = 20; % Fan: fit from data
    end
    
    methods
        function obj = Mccpvd2DofMD()
            obj = obj@Mccpvd2Dof();
            obj.actuator1.Ks = 750;
            obj.actuator2.Ks = 750;
            
            %obj.actuator1.B = 0.03;
            %obj.actuator1.C = 0.125;
            
            %obj.actuator2 = ActMccpvd();
        end
        function xdot = dynamics(model, x, u)
            % uu: (m1 m2 u3)_1 , (m1 m2 u3)_2
            %uu = [x(5:6,:); u(3,:); x(9:10,:); u(6,:)];
            % joint accel
            qddot = model.qddot( x, u);
            % first 4 of xdot
            xdot1 = [x(3:4); qddot];
            % last 8 of xdot
            xdot2 = model.motor_dynamics_2nd(x(5:8,:), u(1:2,:));
            
            xdot3 = model.motor_dynamics_2nd(x(9:12,:), u(4:5,:));
            
            xdot = [xdot1; xdot2; xdot3];
            
        end
        
        function [xdot, xdot_x, xdot_u] = dynamics_with_jacobian_fd(model, x, u)
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
            xdot = model.dynamics(x,u);
            if nargout > 1
                fx = @(x)model.dynamics(x, u);
                dfdx = get_jacobian_fd(fx, x);
                if iscolumn(dfdx), dfdx = dfdx'; end
                % Compute derivative of acceleration w.r.t. u
                %daccdu = zeros(1,model.dimU);
                fu = @(u)model.dynamics(x, u);
                dfdu = get_jacobian_fd(fu, u);
                if iscolumn(dfdu), dfdu = dfdu'; end
                xdot_x = dfdx;
                      
                xdot_u = dfdu;
            end
        end
        
        % mix of analytical and finite diff terms
        function [xdot, xdot_x, xdot_u] = dynamics_with_jacobian_fd_mix(model, x, u)
            % uu: (m1 m2 u3)_1 , (m1 m2 u3)_2
            %uu = [x(5:6,:); u(3,:); x(9:10,:); u(6,:)];
            % joint accel
            qddot = model.qddot( x, u);
            % first 4 of xdot
            xdot1 = [x(3:4); qddot];
            % last 8 of xdot
            [xdot2, xdot2_x, xdot2_u] = model.motor_dynamics_2nd(x(5:8,:), u(1:2,:));
            
            [xdot3, xdot3_x, xdot3_u] = model.motor_dynamics_2nd(x(9:12,:), u(4:5,:));
            
            xdot = [xdot1; xdot2; xdot3];
            
            if nargout > 1
                fx = @(x)model.qddot(x, u);
                daccdx = get_jacobian_fd(fx, x);
                if iscolumn(daccdx), daccdx = daccdx'; end
                % Compute derivative of acceleration w.r.t. u
                %daccdu = zeros(1,model.dimU);
                fu = @(u)model.qddot(x, u);
                daccdu = get_jacobian_fd(fu, u);
                if iscolumn(daccdu), daccdu = daccdu'; end
                xdot_x = [0, 0, 1, 0, zeros(1,8) ;
                          0, 0, 0, 1, zeros(1,8) ;
                          daccdx           ;
                          zeros(4,4), xdot2_x, zeros(4,4);
                          zeros(4,4), zeros(4,4), xdot3_x ];
                      
                xdot_u = [zeros(1,6) ;
                          zeros(1,6) ;
                          daccdu     ;
                          xdot2_u, zeros(4,4);
                          zeros(4,3), xdot3_u, zeros(4,1)];
            end
        end
        
        function acc = qddot(model, x, u)
            uu = [x(5:6,:); u(3,:); x(9:10,:); u(6,:)];
            q = x(1:2,:);
            qdot = x(3:4,:);
            
            tau = model.tau(q, qdot, uu);
            acc = Arm2Dof.compute_qddot(q, qdot, tau, model);
        end
        
        
    
        function [ xdot, xdot_x, xdot_u ] = motor_dynamics_2nd(model, x, u)
            
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
        
        function [ traj ] = reach_ilqr(obj, x0, qf, varargin)
            parser = inputParser();
            
            %addRequired(parser, 'robot');
            addRequired(parser, 'x0');
            addRequired(parser, 'qf');
            addOptional(parser, 'tf', 0.8);
            addOptional(parser, 'dt', 0.01);
            addOptional(parser, 'm2s', [0,0]);
            addOptional(parser, 'sfactor', 1);
            %addOptional(parser, 'ConstraintType', 1);
            addOptional(parser, 'plot', 0);
            addOptional(parser, 'animate', 0);
            addOptional(parser, 'saveanimate', 0);
            addOptional(parser, 'os', 0); % reach in opertional space
            parse(parser, x0, qf, varargin{:});
            
            dt = parser.Results.dt;
            tf = parser.Results.tf;
            m2s = parser.Results.m2s;
            
            t = 0:dt:tf;
            Nt = tf/dt + 1;
            
            
            f = @(x,u)obj.dynamics_with_jacobian_fd( x, u);
            costpara = [];
            
            if parser.Results.os == 0
            
            costpara.Hf = 1e3*diag(ones(1,12));
            %costpara.Hf = 1000*diag([1 1 1 1 0 0 0 0]);
            costpara.Hr = 1*diag([1 1 0.001 1 1 0.001])*parser.Results.sfactor;
            costpara.qf = qf;
            costpara.xf = [ qf(1);qf(2);0;0; qf(1);m2s(1);0;0;qf(2);m2s(2);0;0 ];
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
            opt_param.umax = [6.8 ; 6.8 ; 1 ; 6.8; 6.8; 1] ;
            opt_param.umin = [-6.8; -6.8; 0 ; -6.8;-6.8;0] ;
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
            traj.Jp = (traj.x(:,end)-costpara.xf)'*costpara.Hf*(traj.x(:,end) - costpara.xf);
            else
               traj.y = obj.endpoint(traj.x);
               traj.Jp = (traj.y(:,end)-costpara.yf)'*costpara.Hy*(traj.y(:,end)-costpara.yf);
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
    end
end

