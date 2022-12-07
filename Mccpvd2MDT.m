classdef Mccpvd2MDT < Mccpvd2Dof
    %MCCPVD2MDT Maccepa-vd 2 link. Motor controlled by torque.
    
    properties
        name = 'MACCEPA-VD 2 DoF with 2nd order motor dynamics';
        
        alpha_servo = 20; % Fan: fit from data
    end
    
    methods
        function obj = Mccpvd2MDT()
            %MCCPVD2MDT Construct an instance of this class
            %   Detailed explanation goes here
            obj = obj@Mccpvd2Dof();
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

