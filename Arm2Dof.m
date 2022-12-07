classdef Arm2Dof < Arm
    %ARM2DOF 2 link arm
    %   forward kinematics
    %   x = Arm2Dof.endpoint(q, L)
    %   end point jacobian
    %   J = Arm2Dof.jacobian(q, L)
    % from Arm, inherent step.
    
    properties
        fh
        ph1
        ph2
    end
    
    methods
        function obj = Arm2Dof()
            obj = obj@Arm();
        end
        
        function [] = plot(obj, q)
            [obj.fh, obj.ph1, obj.ph2] = obj.plot_arm(q, obj.L);
            view([-90 90]) 
        end
        
        function [] = update(obj, q)
            obj.update_plot(obj.fh, obj.ph1, obj.ph2, q, obj.L);
        end
        
        function [frames] = animate(obj, x, dt, varargin)
            parser = inputParser();
            
            %addRequired(parser, 'robot');
            addRequired(parser, 'x');
            addRequired(parser, 'dt');
            addOptional(parser, 'goal', [0 0]);
            %addOptional(parser, 'dt', 0.01);
            %addOptional(parser, 'm2f', [0,0]);
            parse(parser, x, dt, varargin{:});
            
            
            
            q = x(1:2,:);
            obj.plot(q(:,1))
            if any( strcmp('goal', varargin))
                goal = parser.Results.goal;
                hold on
                scatter(goal(1),goal(2));
                %hold off
            end
            frames(1) = getframe;
            
            
            for i=2:size(q,2)
                pause(dt)
                obj.update(q(:,i))
                frames(i) = getframe;
                
            end
        end
        
    end
    
    methods (Static)
        
        function qddot = compute_qddot(q, qdot, tau, model)
            Mq = Arm2Dof.compute_Mq(q, model);
            C = Arm2Dof.compute_C(q, qdot, model);
            G = Arm2Dof.compute_G(q, model);
            
            
            
            fc = model.coulomb_friction;
            fv = model.viscous_friction;
            qddot = Mq\(tau - C*qdot - G);
            qddot = qddot - Mq\( fv*qdot + fc*sign(qdot));
            
            %qddot = Mq\(tau - C*qdot - G - fv*qdot - fc*sign(qdot));
        end
        
        % inertia matrix
        function M = compute_Mq(q, model)
            I = model.I; %inertia
            L = model.L; %link length
            M = model.M; % mass
            Lg = model.Lg; % center of mass
            
            M11 = I(1) + I(2) + M(2)*(L(1)^2) + 2*M(2)*L(1)*Lg(2)*cos(q(2));
            M12 = I(2) + M(2)*L(1)*Lg(2)*cos(q(2));
            M21 = I(2) + M(2)*L(1)*Lg(2)*cos(q(2));
            M22 = I(2);
            M = [M11, M12;...
                 M21, M22];
        end
        
        function G = compute_G(q, model)
            L   = model.L;  % link lengths
            M   = model.M;  % mass
            Lg  = model.Lg; % center of gravity
            g   = model.g;

            % gravity vector
            G = -g*[(M(1)*Lg(1)+M(2)*L(1))*cos(q(1))+M(2)*Lg(2)*cos(q(1)+q(2));...
                                        +M(2)*Lg(2)*cos(q(1)+q(2))];
        end
        
        function C = compute_C(q, qdot, model)
            
            L   = model.L;  % link lengths
            M   = model.M;  % mass
            Lg  = model.Lg; % center of gravity

            % Coriolis matrix
            C = M(2)*L(1)*Lg(2)*sin(q(2))*[-2*qdot(2),-qdot(2);...
	                            qdot(1),      0];
        end
        
        %%%% kinematics
        
        % end point position
        function x = endpoint(q, L)
            x1 = L(1)*cos(q(1,:)) + L(2)*cos( q(1,:)+q(2,:) );
            x2 = L(1)*sin(q(1,:)) + L(2)*sin( q(1,:)+q(2,:) );
            
            x = [x1;x2];
        end
        
        % end point jacobian
        function J = jacobian(q, L)
            J(1,1,:) = -L(1)*sin(q(1,:))-L(2)*sin(q(1,:)+q(2,:));
            J(1,2,:) =                  -L(2)*sin(q(1,:)+q(2,:));
            J(2,1,:) =  L(1)*cos(q(1,:))+L(2)*cos(q(1,:)+q(2,:));
            J(2,2,:) =                   L(2)*cos(q(1,:)+q(2,:));
        end
        function dJ = compute_dJ(q,qd, L)
            dJ(1,1,:) = -L(1)*cos(q(1,:))*qd(1,:) - L(2)*cos(q(1,:)+q(2,:))*(qd(1,:) + qd(2,:));
            dJ(1,2,:) =                           - L(2)*cos(q(1,:)+q(2,:))*(qd(1,:) + qd(2,:));
            dJ(2,1,:) = -L(1)*sin(q(1,:))*qd(1,:) - L(2)*sin(q(1,:)+q(2,:))*(qd(1,:) + qd(2,:));
            dJ(2,2,:) =                           - L(2)*sin(q(1,:)+q(2,:))*(qd(1,:) + qd(2,:));
        end
        function ddJ = compute_ddJ(q, qdd, L)
            ddJ(1,1,:) = L(1)*sin(q(1,:))*qdd(1,:) + L(2)*sin(q(1,:)+q(2,:))*(qdd(1,:) + qdd(2,:));
            ddJ(1,2,:) =                           + L(2)*sin(q(1,:)+q(2,:))*(qdd(1,:) + qdd(2,:));
            ddJ(2,1,:) = -L(1)*cos(q(1,:))*qdd(1,:) - L(2)*cos(q(1,:)+q(2,:))*(qdd(1,:) + qdd(2,:));
            ddJ(2,2,:) =                           - L(2)*cos(q(1,:) +q(2,:))*(qdd(1,:) + qdd(2,:));
        end
        function [fh,ph1,ph2] = plot_arm(q, L)
            x11 = L(1)*cos(q(1,:));
            x12 = L(1)*sin(q(1,:));
            x21 = L(1)*cos(q(1,:)) + L(2)*cos( q(1,:)+q(2,:) );
            x22 = L(1)*sin(q(1,:)) + L(2)*sin( q(1,:)+q(2,:) );
            fh = figure;
            LL = L(1)+L(2);
            xlim( [-LL-0.2, LL+0.2] )
            ylim([-LL-0.2,LL+0.2])
            hold on
            %line([0,x11],[0,x12])
            %line([x11,x21],[x12,x22])
            ph1 = plot([0,x11],[0,x12],'-bo','LineWidth',5,'MarkerSize',12,...
                'MarkerFaceColor','k','MarkerEdgeColor','none');
            ph2 = plot([x11,x21],[x12,x22],'-bo','LineWidth',5,'MarkerSize',...
                9,'MarkerFaceColor','k','MarkerEdgeColor','none');
            hold off
        end
        
        function [] = update_plot(h, ph1,ph2, q, L)
            x11 = L(1)*cos(q(1,:));
            x12 = L(1)*sin(q(1,:));
            x21 = L(1)*cos(q(1,:)) + L(2)*cos( q(1,:)+q(2,:) );
            x22 = L(1)*sin(q(1,:)) + L(2)*sin( q(1,:)+q(2,:) );
            
            ph1.XData = [0,x11];
            ph1.YData = [0,x12];
            ph2.XData = [x11,x21];
            ph2.YData = [x12,x22];
        end
        function [] = savegif(frames, gifname)
            for idx = 1:length(frames)
                im = frame2im(frames(idx));
                [A,map] = rgb2ind(im,256);
                if idx == 1
                    imwrite(A,map,gifname,'gif','LoopCount',Inf,'DelayTime',0.02);
                else
                    imwrite(A,map,gifname,'gif','WriteMode','append','DelayTime',0.02);
                end
            end
        end
    end
    
end

