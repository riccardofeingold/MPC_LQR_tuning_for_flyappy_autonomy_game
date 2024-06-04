% LQR Testing
% Q = [
%     1200 0 0 0;
%     0  500 0 0;
%     0 0 1000 0;
%     0 0 0 5
% ];
% 
% R = [
%    	0.001 0;
%     0 0.001
% ];

Q = [
    10 0 0 0;
    0  0.0001 0 0;
    0 0 1000 0;
    0 0 0 15
];

R = [
   	0.001 0;
    0 0.001
];

params = generate_params(Q, R);
Ad = params.model.Ad
Bd = params.model.Bd
Cd = params.model.Cd

N = 30;
mpc = MPC(Q, R, N, params);

SimHorizon = 100;
x0 = [0 0 0.5 0.5]';
x_ref = [4.7, 0.5, 6.0, 0.0]';

x = zeros(4, SimHorizon + 1);
u = zeros(2, SimHorizon - 1);

xk = x0;
x(:, 1) = x0;

[xs, us] = computeSteadyState(x_ref, params)
[uk, objective, feasible] = mpc.eval(xk-xs, xs, us)
for k = 2:SimHorizon+1
    [uk, objective, feasible] = mpc.eval(xk-xs, xs, us);
    xk = Ad * xk + Bd * (uk + us);
    x(:, k) = xk;
    u(:, k-1) = uk + us;
end

figure;

subplot(3, 3, 1);
plot(x(1, :));
title("X");

subplot(3, 3, 2);
plot(x(2, :));
title("X VEL");

subplot(3, 3, 3);
plot(x(3, :));
title("Y");

subplot(3, 3, 4);
plot(x(4, :));
title("Y VEL");

subplot(3, 3, 5);
plot(u(1, :));
title("ax");

subplot(3, 3, 6);
plot(u(2, :));
title("ay");
% MPC Testing
