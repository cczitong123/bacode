classdef Kawato2Dof < Arm2Dof
    %KAWATO2DOF Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        name = 'kawato';
        
        dimQ = 2; % joint Dof
        dimU = 6; % control Dof
        
        % Model geometry and dynamics parameters.
        % Choices for these are based on the paramters in the Katayama/Kawato paper.
        L = 0.3*ones(2,1);
        I = [0.0477;0.0588];
        M = [1.59; 1.44];
        Lg = [0.18; 0.21];
        g = 0;
        
        qmax = pi/2*ones(2,1);
        qmin = -pi/2*ones(2,1); %[-pi*2/3;0]; % what is a good value for this?
        umax = ones(6,1);
        umin = zeros(6,1);
        
        viscous_friction = 0;
        coulomb_friction = 0;
        
        % Muscle parameters vector; note the indices
        % 1 = shoulder flexor
        % 2 = shoulder extensor
        % 3 = elbow flexor
        % 4 = elbow extensor
        % 5 = two-joint flexor
        % 6 = two-joint extensor
        A   = 0.01*[4.0 -4.0 0 0 2.8 -2.8;...  % constant moment arm matrix
                  0 0 2.5 -2.5 3.5 -3.5]';   % negative sign added
                                             % to extensor muscles

        gk  = ones(6,1)*1621.6;                % elasticity coefficients
        k0  = ones(6,1)*810.8;                 % initial elasticity
        gb  = ones(6,1)*108.1;                 % viscosity coefficients
        b0  = ones(6,1)*54.1;                              % initial viscosity 
        gr  = 0.01*[3.491; 3.491; 2.182; 2.182; 5.498; 5.498]; % muscle activation constant

        lm_l0
    end
    
    methods
        function obj = Kawato2Dof()
            obj = obj@Arm2Dof();
            
            % now define lm-l0 according to the standard posture
            A = obj.A;
            K0 = diag(obj.k0);
            P = A'*K0;
            Pinv = P'*inv(P*P'); % pseudo inverse of P=A'*K0
        
            q0 = [45;70]*pi/180;
            delta0 = [0;0;0;0;0;0]; % defines pre-tension
        
            % constant for standard posture
            obj.lm_l0 = Pinv*A'*K0*A*q0 + (eye(6)-Pinv*P)*delta0;
        end
        
        % muscle damping
        function bm = compute_bm(model, u)
            bm = model.b0 + model.gb.*u;
        end
        
        % muscle tension
        function T = compute_Tm(model, q, qdot, u)
            N = size(q,2);
            T = zeros(6,N);
            for n = 1:N
                km = model.compute_km(u(:,n) );
                bm = model.compute_bm(u(:,n) );
                
                ldot = -model.A*qdot(:,n);
                
                T(:,n) = km.*( model.lm_l0 - model.A*q(:,n) + model.gr.*u(:,n) ) + bm.*ldot;
            end
        end
        
        function km = compute_km(model, u)
            km = model.k0 + model.gk.*u;
            
        end
        
        function tau = tau(model, q, qdot, u)
            T = model.compute_Tm(q, qdot, u);
            tau = model.A'*T;
        end
        
        function qddot = qddot(model, q, qdot, u)
            tau = model.tau(q, qdot, u);
            qddot = Arm2Dof.compute_qddot(q, qdot, tau, model);
        end
        
        function xdot = dynamics(model, x, u)
            qddot = model.qddot( x(1:2), x(3:4), u);
            xdot = [x(3:4); qddot];
        end
        
        function x = endpoint(model, q)
            x = Arm2Dof.endpoint(q, model.L); 
        end
        
        function J = jacobian(model, q)
            J = Arm2Dof.jacobian(q, model.L);
        end
    end
    
end

