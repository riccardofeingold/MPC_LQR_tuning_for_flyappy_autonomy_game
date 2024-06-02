function [params] = generate_params()
    % Define System Dynamics
    params = struct();
    Ts = 1/30;
    
    params.model = struct(...
        'nx', 4, ...
        'nu', 2, ...
        'ny', 2, ...
        'Ts', Ts ...
        );
   params.model.A = [
        0 1 0 0;
        0 0 0 0;
        0 0 0 1;
        0 0 0 0
    ];
   params.model.B = [
        0 0;
        1 0;
        0 0;
        0 1
    ];
   params.model.C =  [
        1 0 0 0;
        0 0 1 0
   ];

    sys = ss(params.model.A, params.model.B, eye(size(params.model.A, 1)), zeros(size(params.model.A, 1), size(params.model.B, 2)));
    sysd = c2d(sys, Ts);

    params.model.Ad = sysd.A;
    params.model.Bd = sysd.B;
    params.model.Cd = params.model.C;

    % Define constraints
    params.constraints.Hx = [
        0 0 1 0;
        0 0 -1 0
    ];

    yUpperBound = 2;
    yLowerBound = -2;
    params.constraints.hx = [
        yUpperBound;
        -yLowerBound
    ];

    params.constraints.Hu = [
        1 0;
        -1 0;
        0 1;
        0 -1
    ];

    axUpperBound = 3;
    axLowerBound = -axUpperBound;
    ayUpperBound = 35;
    ayLowerBound = -ayUpperBound;
    params.constraints.hu = [
        axUpperBound;
        -axLowerBound;
        ayUpperBound;
        -ayLowerBound
    ];
end