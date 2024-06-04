function [x_s, u_s] = computeSteadyState(xref, params)
    C = [
        1 0 0 0;
        0 0 1 0
    ];
    C = eye(4);
    state_space = [
        eye(size(params.model.Ad, 1)) - params.model.Ad, params.model.Bd;
        C, zeros(size(C, 1), size(params.model.B, 2))
    ];
    
    RHS = [
        zeros(params.model.nx, 1);
        xref
    ];

    result = state_space \ RHS;

    x_s = result(1:params.model.nx);
    u_s = result(params.model.nx+1:end);
end