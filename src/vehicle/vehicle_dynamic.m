function dx = vehicle_dynamic(x, u, vp, prm)
% u = [F_drive; delta]
[ddx, ddy, ddphi] = vehicle_model(x(4), x(5), u(1), u(2), ...
    vp.h_a, vp.h_b, vp.m, vp.track, x(3), x(6), prm);
dx = [x(4); x(5); x(6); ddx; ddy; ddphi];
end