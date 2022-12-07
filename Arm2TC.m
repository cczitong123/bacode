classdef Arm2TC < Arm2Dof
    % ARM2TC torque controlled 2-link arm
    % state: q1 q2 q1dot q2dot
    properties (Access = private)
      name = 'Arm2 torque controlled'
    end
    properties
        
        
        dimQ = 2 % 2 joints
        dimU = 2 % 2 torques
        
        L = [0.3;0.33];
        I = [0.025;0.045]; 
        M = [1.4; 1]; %[1.59; 1.44];
        Lg = [0.11; 0.16];
        g = 0;
        
        viscous_friction = 0.05;
        coulomb_friction = 0;
        
        %%%%
        
        
        %%%%
    end
    
    methods
        function obj = Arm2TC()
            %ARM2TC Construct an instance of this class
            %   Detailed explanation goes here
            obj = obj@Arm2Dof();
        end
        
        function qddot = qddot(obj, q, qdot, u)
            % u: (m1 m2 u3)_1 , (m1 m2 u3)_2
            
            qddot = Arm2Dof.compute_qddot(q, qdot, u, obj);
        end
        
        % x: q1 q2 qd1 qd2
        % u: u
        function xdot = dynamics(obj, x, u)
            qddot = obj.qddot( x(1:2,:), x(3:4,:), u);
            xdot = [x(3:4); qddot];
        end
        
        function [xdot, xdot_x, xdot_u] = dynamics_with_jacobian(obj,x,u)
            xdot = obj.dynamics(x, u);
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
        
        function x = endpoint(obj, q)
            x = Arm2Dof.endpoint(q, obj.L); 
        end
        
        function J = jacobian(obj, q)
            J = Arm2Dof.jacobian(q, obj.L);
        end
        
        function [l, l_x, l_xx, l_u, l_uu, l_ux] = costfnt_reach(obj, x,u,t, costpara)
            %Rtf  = diag([w_tf, 0, 0, 0, 0, 0]);
            Hf = costpara.Hf;
            Hr = costpara.Hr;
            xf = costpara.xf;
            if (isnan(u))
                % final cost
                %fh = @(x)obj.costf_reach(x);
                
                l = ((x-xf)')*Hf*(x-xf)/2;
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
                    
                    l_x = [0 0 0 0];
                    l_u = u'*Hr;
                    l_xx = zeros(4,4);
                    l_ux = zeros(2,4);
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
            %addOptional(parser, 'm2min', [0,0]);
            %addOptional(parser, 'sfactor', 3);
            %addOptional(parser, 'ConstraintType', 1);
            addOptional(parser, 'plot', 0);
            addOptional(parser, 'animate', 0);
            parse(parser, x0, qf, varargin{:});
            
            dt = parser.Results.dt;
            tf = parser.Results.tf;
            t = 0:dt:tf;
            Nt = tf/dt + 1;
            
            
            f = @(x,u)obj.dynamics_with_jacobian( x, u);
            
            costpara = [];
            costpara.Hf = 10000*diag([1 1 1 1]);
            costpara.Hr = 1000*diag([1 1]);
            costpara.qf = qf;
            costpara.xf = [ qf(1);qf(2);0;0 ];
            costpara.tf = tf;
            
            %costpara.L = robot.L ; % link lengths
            %%%%
            j = @(x,u,t)obj.costfnt_reach(x, u, t, costpara);
            %%%% Define parameters of optimisation 
            opt_param = [];
            
            opt_param.lambda_init = 1;
            opt_param.lambda_max  = 1e10;
            opt_param.iter_max = 100;
            opt_param.online_plotting = 0;
            opt_param.online_printing = 2;
            opt_param.dcost_converge = 1e-3;
            opt_param.solver = 'rk4';
            
            opt_param.T = tf;
            %%%
            opt_param.umax = [100 ; 100 ] ;
            opt_param.umin = [-100; -100] ;
            % u0 can be full command sequence or just initial point
            %u0 = [qf(1); m2min(1)+0.1; 0; qf(2) ;m2min(2)+0.1 ; 0];
            
            u0 = rand(2,1);
            %%%%
            tic
            %traj = ILQRController.ilqr(f, j, dt, Nt, x0, u0, opt_param);
            
            traj = iLQRv1(f, j, dt, Nt, x0, u0, opt_param);
            traj.t = t;
            %traj = val_traj_mccpvd(robot, task, traj);
            %if nargin == 0
            %    plot_traj_mccpvd1(traj);
            %end
            
            toc
            disp('ILQR done')
            
            if parser.Results.plot == 1
                figure
                subplot(3,1,1)
                plot(traj.t, traj.x(1,:))
                title('Joint 1')
                ylabel('q_1 [rad]')
                subplot(3,1,2)
                plot(traj.t, traj.x(2,:))
                title('Joint 2')
                ylabel('q_2 [rad]')
                xlabel('t [s]')
                subplot(3,1,3)
                plot(traj.t(1:end-1), traj.u )
            end
            
            if parser.Results.animate == 1
                frames = obj.animate(traj.x,0.01);
                Arm2Dof.savegif(frames, '2link_reaching.gif')
            end
        end
        
        function C = costfun_reachFast()
            
        end
    end
end

