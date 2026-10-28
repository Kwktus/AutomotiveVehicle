function vp = vehicle_params()
vp.m     = 1500;     % масса, кг
vp.h_a   = 1.2;      % ЦМ -> передняя ось, м
vp.h_b   = 1.4;      % ЦМ -> задняя ось, м
vp.track = 1.6;      % колея, м
vp.v0    = 0;       % начальная скорость, м/с

vp.prm = struct('mu',0.8, 'Cy',45000, 'h_g',0.5, 'f_roll',0.015);
end