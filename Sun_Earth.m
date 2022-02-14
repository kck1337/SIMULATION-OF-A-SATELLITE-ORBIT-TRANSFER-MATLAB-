close all;
orbit = 1.3e3;
load traj.mat
steps = 1000 ;
k=1;
j=1;
ang = linspace(0,360,(length(States_Y_E)/steps)+1);
for i = 1:steps:length(States_Y_E)
%      f = figure();
    if j==length(States_Y_S)
        j = 1
    end
    clf
    sun()
    hold on
    plot3(States_X_E,States_Y_E,States_Z_E)
    view([1 -1 1])
    axis equal
    earth(i,ang(k))
%     moon(i,j,ang(k),orbit)
   circle(States_X_E(i),States_Y_E(i),States_Z_E(i),orbit,ang(k),j);
    hold off
    refreshdata
    drawnow
    k=k+1;
    j = j+1;
%      saveas(f,sprintf("  FIG%d.png",k));
     %set (position(i-1), 'Visible', 'off')   
%     xlim([-inf inf])
%     ylim([-inf inf])
end
%hold off


function ot = circle(x,y,z,r,ang,j)
%x and y are the coordinates of the center of the circle
%r is the radius of the circle
%0.01 is the angle step, bigger values will draw the circle faster but
%you might notice imperfections (not very smooth)
load traj.mat
xp=States_X_S*r;
yp=States_Y_S*r;
zp=States_Z_S*r;

ot = plot3(xp+x,yp+y,zp+z,'r-');
rotate(ot,[0 0 1],2*ang,[x y z])
end

function moon(i,j,ang,r)
load traj.mat
[m1,m2,m3] = sphere(40);
r_moon= 1737.44;
m1=m1*r_moon; m2 = m2*r_moon; m3 = m3*r_moon;
moon = surface(m1+States_X_E(i)+States_X_S(j)*r,m2+States_Y_E(i)+States_Y_S(j)*r,m3+States_Z_E(i)+States_Z_S(j)*r);
rotate(moon,[0 0 1],2*ang,[States_X_E(i) States_Y_E(i) States_Z_E(i)])
moon.FaceColor = 'black';
axis equal
end

function earth(i,ang)
[e1,e2,e3] = sphere;
load topo
load traj.mat
r_earth= 1e3*6371;
e1=e1*r_earth; e2 = e2*r_earth; e3 = e3*r_earth;
earth = surface(e1+States_X_E(i),e2+States_Y_E(i),e3+States_Z_E(i));
% earth = surface(e1,e2,e3);
rotate(earth,[0 0.3 1],6*ang,[States_X_E(i) States_Y_E(i) States_Z_E(i)])
earth.FaceColor = 'texturemap';
earth.CData = topo;
earth.EdgeColor = 'none';
earth.FaceLighting = 'gouraud';
earth.SpecularStrength = 0.4;
axis equal
end

function sun()
[s1,s2,s3] = sphere(50);
r = 1.2e1*695508;
sun = surface(s1*r,s2*r,s3*r);
sun.EdgeColor = 'none';
map = imread('Images/Sun_map.jpg'); % read planet map
map = imrotate(map,180); % rotate upright
Cord = [0 0 0 r];
warp(s1*Cord(1,4),s2*Cord(1,4),s3*Cord(1,4),map);
view ([20 20])
axis equal
end

