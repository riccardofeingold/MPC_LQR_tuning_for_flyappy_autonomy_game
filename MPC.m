classdef MPC
    properties
        yalmip_optimizer
    end

    methods
        function obj = MPC(Q, R, N, params)
            A = params.model.Ad;
            B = params.model.Bd;
            Hu = params.constraints.Hu;
            hu = params.constraints.hu;
            Hx = params.constraints.Hx;
            hx = params.constraints.hx;

            U = sdpvar(repmat(params.model.nu, 1, N), repmat(1, 1, N));
            X = sdpvar(repmat(params.model.nx, 1, N+1), repmat(1, 1, N+1));
            X0 = sdpvar(params.model.nx, 1);

            constraints = [X{1} == X0, X{N+1} == zeros(params.model.nx, 1)];
            objective = 0;

            for k = 1:N
                objective = objective + X{k}' * Q * X{k} + U{k}' * R * U{k};
                % constraints = [constraints, X{k}(1, 1) <= X{k+1}(1, 1), X{k+1} == A*X{k} + B*U{k}, Hx * X{k+1} <= hx, Hu * U{k} <= hu];
                constraints = [constraints, X{k+1} == A*X{k} + B*U{k}, Hu * U{k} <= hu];
            end

            [~, P_inf, ~] = dlqr(A, B, Q, R);
            objective = objective + X{N+1}' * P_inf * X{N+1};
            
            opts = sdpsettings('verbose',1,'solver','quadprog');
            obj.yalmip_optimizer = optimizer(constraints,objective,opts,X0,{U{1} objective});
        end

        function [u, objective, feasible] = eval(obj, x)
            [optimizer_out,errorcode] = obj.yalmip_optimizer(x);

            [u, objective] = optimizer_out{:};

            feasible = true;
            if (errorcode ~= 0)
                feasible = false;
            end
        end
    end
end