function truth = run_truth(vp, input_fun, t_end, dt)
N = floor(t_end/dt) + 1;
t = (0:N-1)'*dt;

x       = [0; 0; 0; vp.v0; 0; 0];     % стартовое состояние
ay_prev = 0;
prm     = vp.prm;

S = zeros(N,6);  A = zeros(N,4);      % состояния и ускорения

for k = 1:N
    u = input_fun(t(k), x);
    prm.a_y_prev = ay_prev;           % заморожено на шаг
    f = @(tt, xx) vehicle_dynamic(xx, u, vp, prm);

    dx = f(t(k), x);                  % ускорения в текущем состоянии
    phi = x(3);
    ax_b =  dx(4)*cos(phi) + dx(5)*sin(phi);
    ay_b = -dx(4)*sin(phi) + dx(5)*cos(phi);

    S(k,:) = x';
    A(k,:) = [dx(4) dx(5) ax_b ay_b];
    ay_prev = ay_b;

    if k < N
        x = rk4_step(f, t(k), x, dt);
    end
end

truth.t = t;
truth.x = S(:,1);  truth.y = S(:,2);  truth.phi = S(:,3);
truth.vx = S(:,4); truth.vy = S(:,5); truth.dphi = S(:,6);
truth.ax_w = A(:,1); truth.ay_w = A(:,2);
truth.ax_b = A(:,3); truth.ay_b = A(:,4);
end