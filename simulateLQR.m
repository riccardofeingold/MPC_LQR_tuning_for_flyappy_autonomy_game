params = generate_params();
Ad = params.model.Ad
Bd = params.model.Bd
Cd = params.model.Cd

% LQR Testing
Q = [
    100 0 0 0;
    0  1 0 0;
    0 0 300 0;
    0 0 0 1
];

R = [
   	50 0;
    0 1
];

lqr_controller = LQR(Q, R, Ad, Bd);

K = lqr_controller.K

SimHorizon = 100;
x0 = [1 1 0 0]';
x_ref = [2, 0, -2, 0]';
x = zeros(4, SimHorizon + 1);
u = zeros(2, SimHorizon - 1);

xk = x0;
x(:, 1) = x0;

for k = 2:SimHorizon+1
    uk = lqr_controller.eval(xk - x_ref);
    xk = Ad * xk + Bd * uk;
    x(:, k) = xk;
    u(:, k-1) = uk;
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
