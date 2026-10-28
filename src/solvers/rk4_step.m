function x_next = rk4_step(f, t, x, dt)
% f  - handle производной: dx = f(t, x)
k1 = f(t,        x);
k2 = f(t + dt/2, x + dt/2*k1);
k3 = f(t + dt/2, x + dt/2*k2);
k4 = f(t + dt,   x + dt*k3);
x_next = x + dt/6*(k1 + 2*k2 + 2*k3 + k4);
end