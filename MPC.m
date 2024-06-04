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
            H = params.constraints.H;
            h = params.constraints.h;
            
            nx = params.model.nx;
            nu = params.model.nu;

            U = sdpvar(repmat(params.model.nu, 1, N), repmat(1, 1, N));
            X = sdpvar(repmat(params.model.nx, 1, N+1), repmat(1, 1, N+1));
            XS = sdpvar(nx, 1);
            US = sdpvar(nu, 1);
            X0 = sdpvar(params.model.nx, 1);

            % get delta formulation
            delta_hx = hx - Hx * XS;
            delta_hu = hu - Hu * US;
            
            Hf = [
                % 1, 0, 0, 0;
                % -1, 0, 0, 0;
                % 0, 1, 0, 0;
                % 0, -1, 0, 0;
                0, 0, 1, 0;
                0, 0, -1, 0;
                0, 0, 0, 1;
                0, 0, 0, -1
            ];

            hf = zeros(size(Hf, 1), 1);
            
            objective = 0;
            constraints = [X{1} == X0, Hf * X{N+1} <= hf];

            for k = 1:N
                objective = objective + X{k}' * Q * X{k} + U{k}' * R * U{k};
                constraints = [constraints, X{k+1} == A*X{k} + B*U{k}, Hx * X{k+1} <= delta_hx, Hu * U{k} <= delta_hu];
            end
            
            % constraints = [X{1} == X0];
            % 
            % for k = 1:N
            %     objective = objective + X{k}' * Q * X{k} + U{k}' * R * U{k};
            %     constraints = [constraints, X{k+1} == A*X{k} + B*U{k}, Hu * U{k} <= delta_hu];
            % end

            [~, P_inf, ~] = dlqr(A, B, Q, R);
            objective = objective + X{N+1}' * P_inf * X{N+1};
            
            opts = sdpsettings('verbose',1,'solver','quadprog');
            obj.yalmip_optimizer = optimizer(constraints,objective,opts,{X0, XS, US},{U{1} objective});
        end

        function [u, objective, feasible] = eval(obj, x, xs, us)
            [optimizer_out,errorcode] = obj.yalmip_optimizer(x, xs, us);

            [u, objective] = optimizer_out{:};

            feasible = true;
            if (errorcode ~= 0)
                feasible = false;
            end
        end
    end
end