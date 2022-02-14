mu = 398600.436233; % Gravitational Parameter for earth
req = 6378.1363;    % Equitorial Radius of Earth
global hn1 hn2 hn3  % normalized radii
% user inputs
while (1)
       fprintf('\n\n Input the initial altitude (kilometers)\n');
        alt1 = input('? ');
        if (alt1 > 0.0)
        break;
    end
end
while (1)
       fprintf('\n\n Input the final altitude (kilometers)\n');
       alt2 = input('? ');
        if (alt2 > 0.0)
        break;
    end
  end
% Radius of initial and final orbits from the center of earth in KM
r1 = req + alt1;
r2 = req + alt2;
% circular velocity of the initial and final orbits km/sec
v1 = sqrt(mu / r1);
v2 = sqrt(mu / r2);
% Semi major axis of transfer orbit in km
smat = 0.5 * (r1 + r2);
% eccentricity of transfer orbit
e = (max(r1, r2) - min(r1, r2)) / (r1 + r2);
% transfer orbit perigee and apogee radii and velocities
rp = smat * (1.0 - e);
ra = smat * (1.0 + e);
vt1 = sqrt(2.0 * mu * ra / (rp * (rp + ra)));
vt2 = sqrt(2.0 * mu * rp / (ra * (rp + ra)));
% orbit transfer
if (r2 > r1)
    % higher-to-lower transfer
    dv1 = vt1 - v1;
    dv2 = v2 - vt2;
else
    % lower-to-higher transfer
    dv1 = v1 - vt2;
    dv2 = vt1 - v2;
end
% load orbital elements arrays, create state vectors and plot orbits
oevi(1) = r1;
oevi(2) = 0.0;
oevi(3) = 0.0;
oevi(4) = 0.0;
oevi(5) = 0.0;
% Finding true anomaly (radians)
if (alt2 > alt1)
        oevi(6) = 0.0;
    else
        oevi(6) = deg2rad(180.0);
    end
[ri, vi] = orb2eci(mu, oevi);
oevti(1) = smat;
oevti(2) = e;
oevti(3) = 0.0;
oevti(4) = 0.0;
oevti(5) = 0.0;
% determine correct true anomaly (radians)
if (alt2 > alt1)
        oevti(6) = 0.0;
    else
        oevti(6) = deg2rad(180.0);
    end
[rti, vti] = orb2eci(mu, oevti);
oevtf(1) = smat;
oevtf(2) = e;
oevtf(3) = 0.0;
oevtf(4) = 0.0;
oevtf(5) = 0.0;
% determine correct true anomaly (radians)
if (alt2 > alt1)
        oevtf(6) = deg2rad(180.0);
    else
        oevtf(6) = 0.0;
    end
[rtf, vtf] = orb2eci(mu, oevtf);
oevf(1) = r2;
oevf(2) = 0.0;
oevf(3) = 0.0;
oevf(4) = 0.0;
oevf(5) = 0.0;
% determine correct true anomaly (radians)
if (alt2 > alt1)
        oevf(6) = deg2rad(180.0);
    else
        oevf(6) = 0.0;
end
[rf, vf] = orb2eci(mu, oevf);
% compute orbital periods
tp1 = 2.0 * pi * oevi(1) * sqrt(oevi(1) / mu);
tp2 = 2.0 * pi * oevti(1) * sqrt(oevti(1) / mu);
tp3 = 2.0 * pi * oevf(1) * sqrt(oevf(1) / mu);
dt1 = tp1 / 360;
t1 = -dt1;
dt2 = 0.5 * tp2 / 360;
t2 = -dt2;
dt3 = tp3 / 360;
t3 = -dt3;
for i = 1:1:361
        t1 = t1 + dt1;
        t2 = t2 + dt2;
        t3 = t3 + dt3;
        % using 2 body initial value problem
        % compute initial orbit "normalized" position vector
        [rwrk, vwrk] = twobody2 (mu, t1, ri, vi);
        rp1_x(i) = rwrk(1) / req;
        rp1_y(i) = rwrk(2) / req;
        rp1_z(i) = rwrk(3) / req;
        % compute transfer orbit position vector
        [rwrk, vwrk] = twobody2 (mu, t2, rti, vti);
        rp2_x(i) = rwrk(1) / req;
        rp2_y(i) = rwrk(2) / req;
        rp2_z(i) = rwrk(3) / req;
        % compute final orbit position vector
        [rwrk, vwrk] = twobody2 (mu, t3, rf, vf);
        rp3_x(i) = rwrk(1) / req;
        rp3_y(i) = rwrk(2) / req;
        rp3_z(i) = rwrk(3) / req;
    end
hold on;
% plot earth
earth_map = imread('Images/Earth_map.jpg'); % read planet map
earth_map = imrotate(earth_map,180);        % rotate upright
[x y z] = sphere(24);
load topo
h = surf(x, y, z);
h.FaceColor = 'texturemap';
h.CData = topo;
h.EdgeColor = 'none';
h.FaceLighting = 'gouraud';
h.SpecularStrength = 0.4;
% plot initial orbit
plot3(rp1_x, rp1_y, rp1_z, '-r', 'LineWidth', 1.5);
plot3(rp1_x(1), rp1_y(1), rp1_z(1), 'ob');
% plot transfer orbit
plot3(rp2_x, rp2_y, rp2_z, '-b', 'LineWidth', 1.5);
plot3(rp2_x(end), rp2_y(end), rp2_z(end), 'ob');
% plot final orbit
plot3(rp3_x, rp3_y, rp3_z, '-g', 'LineWidth', 1.5);
xlabel('X coordinate', 'FontSize', 12);
ylabel('Y coordinate', 'FontSize', 12);
zlabel('Z coordinate', 'FontSize', 12);
axis equal;
view(50, 20);
rotate3d on;