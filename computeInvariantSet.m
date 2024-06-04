function [H, h] = computeInvariantSet(Q, R, A, B, Hx, hx, Hu, hu)
    lqr_controller = LQR(Q, R, A, B);
    lqr_system = LTISystem('A', A - B*lqr_controller.K);
    Xp = Polyhedron('A', [Hx; -Hu * lqr_controller.K], 'b', [hx; hu]);
    lqr_system.x.with('setConstraint');
    lqr_system.x.setConstraint = Xp;

    inv_set_lqr = lqr_system.invariantSet();
    H = inv_set_lqr.A;
    h = inv_set_lqr.b;
end