classdef Arm < handle
    %ARM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function model = Arm()
            
        end
        function [xn] = step(model, x, u, dt, sdt)
            %tic
            %sdt = 0.001;
            t = 0;
            
            while t < dt
            g1 = sdt*model.dynamics(x            ,u);
            g2 = sdt*model.dynamics(x+.5*g1,u);
            g3 = sdt*model.dynamics(x+.5*g2,u);
            g4 = sdt*model.dynamics(x+   g3,u);
            x = x + (1/6)*(g1 + 2*g2 + 2*g3 + g4);
            t = t + sdt;
            end
            xn = x;
            %toc
        end
        
        function x = simulate_feedforward ( model, varargin )
            x0 = varargin{1}; u = varargin{2}; dt = varargin{3};
            if length(varargin) == 4
                sdt = varargin{4};
            else
                sdt = 0.005;
            end
            
            % dt: control timestep
            N = size(u,2)+1; % N of time index
            x = nan(size(x0,1),N); x(:,1)=x0; % initialise x
            for n=1:N-1
                x(:,n+1) = model.step ( x(:,n), u(:,n), dt, sdt );
            end

        end
    end
    
end

