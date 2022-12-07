classdef Qbmove
    %QBMOVE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        k1 =  0.0227 %0.0026 % Nm 
        k2 =  0.0227 %0.0216 %0.0011 % Nm 
        a1 =  6.7328 %8.9995 %
        a2 =  6.7328 %6.9602 %8.9989 %
    end
    
    methods
        function obj = Qbmove(varargin)
        
        end
        
        function tau = torque(obj, x, q1, q2)
            tau = obj.k1*sinh(obj.a1*(q1 - x)) + obj.k2*sinh(obj.a2*(q2-x));
        end
        
        function s = stiffness(obj, x, q1, q2)
            s = obj.a1*obj.k1*cosh(obj.a1*(x-q1)) + obj.a2*obj.k2*cosh(obj.a2*(x-q2));
        end
        
        function E = energy(obj, x, q1, q2)
            E = obj.k1*(cosh(obj.a1*(x-q1))-1)/obj.a1 + ...
                obj.k2*(cosh(obj.a2*(x-q2))-1)/obj.a2;    
        end
        
        function ep = ep(~, q1, q2)
            ep = (q1 + q2)/2;
        end
    end
    
end

