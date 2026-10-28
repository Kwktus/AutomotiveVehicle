function [ddx, ddy, dd_phi] = vehicle_model(v_x, v_y, F_drive, delta, ...
        h_a, h_b, m, track, phi, dphi, prm)
% v_x,v_y  - скорость ЦМ в МИРОВОЙ СК, м/с
% F_drive  - суммарная тяга задней оси, Н (или тормозное усилие со знаком)
% delta    - угол поворота передних колёс, рад
% h_a,h_b  - ЦМ -> передняя/задняя ось, м
% track    - КОЛЕЯ, м   (не путать с базой L_wb = h_a+h_b)
if nargin < 11 || isempty(prm), prm = struct(); end
g = getf(prm,'g',9.81);  h_g = getf(prm,'h_g',0.5);
mu = getf(prm,'mu',0.8); f_roll = getf(prm,'f_roll',0.015);
Cy = getf(prm,'Cy',45000);          % Н/рад на колесо
rho= getf(prm,'rho',1.225); Cd=getf(prm,'Cd',0.32); A=getf(prm,'A',2.2);
k_roll_f = getf(prm,'k_roll_f',0.5); % доля поперечного переноса на переднюю ось
a_y_prev = getf(prm,'a_y_prev',0);
v_eps = getf(prm,'v_eps',0.5);

L_wb = h_a + h_b;
J    = m/12*(L_wb^2 + track^2);      % ИСПРАВЛЕНО: реальные габариты

% --- тело ---
Vxb =  v_x*cos(phi) + v_y*sin(phi);
Vyb = -v_x*sin(phi) + v_y*cos(phi);

% --- предварительная продольная перегрузка (для переноса нагрузки) ---
F_drag = 0.5*rho*Cd*A*Vxb*abs(Vxb);
a_x0   = (F_drive - F_drag - f_roll*m*g) / m;

Fz_f = max(0, m*g*h_b/L_wb - m*a_x0*h_g/L_wb);   % защита от отрицательных
Fz_r = max(0, m*g*h_a/L_wb + m*a_x0*h_g/L_wb);

% --- поперечное перераспределение (раньше отсутствовало) ---
dFz   = m*a_y_prev*h_g/track;
Fz_FL = max(0, 0.5*Fz_f - 0.5*dFz*k_roll_f);
Fz_FR = max(0, 0.5*Fz_f + 0.5*dFz*k_roll_f);
Fz_RL = max(0, 0.5*Fz_r - 0.5*dFz*(1-k_roll_f));
Fz_RR = max(0, 0.5*Fz_r + 0.5*dFz*(1-k_roll_f));

% --- скорости точек контакта (x вперёд, y влево) ---
Vx_FL = Vxb - dphi*track/2;  Vy_FL = Vyb + dphi*h_a;
Vx_FR = Vxb + dphi*track/2;  Vy_FR = Vyb + dphi*h_a;
Vx_RL = Vxb - dphi*track/2;  Vy_RL = Vyb - dphi*h_b;
Vx_RR = Vxb + dphi*track/2;  Vy_RR = Vyb - dphi*h_b;

[Fx_FL,Fy_FL] = tire(0,          Vx_FL,Vy_FL, delta, Cy,mu,Fz_FL,f_roll,v_eps);
[Fx_FR,Fy_FR] = tire(0,          Vx_FR,Vy_FR, delta, Cy,mu,Fz_FR,f_roll,v_eps);
[Fx_RL,Fy_RL] = tire(F_drive/2,  Vx_RL,Vy_RL, 0,     Cy,mu,Fz_RL,f_roll,v_eps);
[Fx_RR,Fy_RR] = tire(F_drive/2,  Vx_RR,Vy_RR, 0,     Cy,mu,Fz_RR,f_roll,v_eps);

% --- поворот шинных сил в тело ---
c=cos(delta); s=sin(delta);
Fx_FLb = Fx_FL*c - Fy_FL*s;  Fy_FLb = Fx_FL*s + Fy_FL*c;
Fx_FRb = Fx_FR*c - Fy_FR*s;  Fy_FRb = Fx_FR*s + Fy_FR*c;
Fx_RLb = Fx_RL;              Fy_RLb = Fy_RL;
Fx_RRb = Fx_RR;              Fy_RRb = Fy_RR;

% --- моменты M_z = x*Fy - y*Fx ---
M =  ( h_a*Fy_FLb - ( track/2)*Fx_FLb) + ( h_a*Fy_FRb - (-track/2)*Fx_FRb) ...
   + (-h_b*Fy_RLb - ( track/2)*Fx_RLb) + (-h_b*Fy_RRb - (-track/2)*Fx_RRb);

Fbx = Fx_FLb+Fx_FRb+Fx_RLb+Fx_RRb - F_drag;
Fby = Fy_FLb+Fy_FRb+Fy_RLb+Fy_RRb;

Fwx = Fbx*cos(phi) - Fby*sin(phi);
Fwy = Fbx*sin(phi) + Fby*cos(phi);

ddx = Fwx/m; ddy = Fwy/m; dd_phi = M/J;
end

% ---------------- шина: гладкое насыщение, без разрыва ----------------
function [Fx,Fy] = tire(Fx_req, Vx, Vy, delta, Cy, mu, Fz, f_roll, v_eps)
Fz  = max(Fz, 1e-3);              % защита от деления на 0
Fmax= max(mu*Fz, 1e-3);

alpha = delta - atan2(Vy, abs(Vx) + v_eps);
alpha = max(min(alpha, 1.2), -1.2);      % ограничение tan(alpha)

Fx0 = Fx_req - f_roll*Fz*tanh(Vx/0.5);
Fy0 = Cy*tan(alpha);

% гладкое круговое ограничение (C1, вместо if/else)
n  = 4;
sq = (Fx0/Fmax)^2 + (Fy0/Fmax)^2;
sc = (1 + sq^(n/2))^(-1/n);
Fx = Fx0*sc; Fy = Fy0*sc;
end

function v = getf(s,name,def)
if isfield(s,name), v = s.(name); else, v = def; end
end
