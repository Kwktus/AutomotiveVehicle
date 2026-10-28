close all;
rng(1);
vp = vehicle_params();

v  = vp.v0;                       % 10 м/с
R  = 50;                          % желаемый радиус, м
L  = vp.h_a + vp.h_b;
Cy = vp.prm.Cy;
Wf = vp.m*9.81*vp.h_b/L;
Wr = vp.m*9.81*vp.h_a/L;
Kus = Wf/(2*Cy) - Wr/(2*Cy);      % недостаточная поворачиваемость, рад/g
delta0 = L/R + Kus*v^2/(R*9.81);  % ≈ 0.055 рад
F0 = 700;                         % Н

input_fun = @(t, x) [F0; delta0];
truth = run_truth(vp, input_fun, 60, 0.001);   % полный круг ≈ 31.4 с

figure;
subplot(2,2,1); plot(truth.x, truth.y); axis equal; grid on
xlabel('x, м'); ylabel('y, м'); title('Траектория');
subplot(2,2,2); plot(truth.t, truth.vx, truth.t, truth.vy); grid on
legend('v_x','v_y'); title('Скорости (мировая СК)');
subplot(2,2,3); plot(truth.t, truth.dphi); grid on
title('Угловая скорость, рад/с');
