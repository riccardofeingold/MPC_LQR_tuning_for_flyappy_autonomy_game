classdef LQR
    properties
        K
    end
    
    methods
        function obj = LQR(Q, R, Ad, Bd)
            [obj.K, ~, ~] = dlqr(Ad, Bd, Q, R);
        end
        
        function u = eval(obj, x)
            u = -obj.K * x;
        end
    end
end