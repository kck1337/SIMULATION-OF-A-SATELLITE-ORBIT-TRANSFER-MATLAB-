clear; close all; clc;

%% Select orbit type
Orbit_type = 'Moon_orbit' % types of orbits are Moon_orbit, Low Earth Orbit, High Earth Orbit, Molniya Orbit, Polar Orbit

%% Creating an object 'Sun' with all the relevant properties to calculate the Gravitational effect of the Sun
u_sun = 1.327124400189e11;              % km^3/s^2 Gravitational Constant 
sun_map = imread('Images/Sun_map.jpg'); % read planet map
sun_map = imrotate(sun_map,180);        % rotate upright
sun_Radius = 695508;                    % km
sun_Mass = 1.989e30;                    % kg
Dist_to_Sun_self = 1;                        % km
%J_num = [0];

%% Creating an object 'Earth' with all the relevant properties to calculate the orbital trajectories
u_earth = 3.986004415e5;                    % km^3/s^2
earth_map = imread('Images/Earth_map.jpg'); % read planet map
earth_map = imrotate(earth_map,180);        % rotate upright
earth_Radius = 6371;                        % km
earth_Mass = 5.9724e24;                    % kg
Dist_to_Sun = 149.6e6;                      % km
%J_num = [0;0.0010826269;-0.0000025323;-0.0000016204]; % zonal harmonic coefficients 

%% Creating an object 'Moon' for visualization of the body as a satellite
moon_map = imread('Images/Moon_map.jpg');
moon_map = imrotate(moon_map,180);
moon_radius = 1738.1;            %km
moon_mass = 7346e19;            %kg
% Dist_to_Earth = 3.78e5;           %km

%% Taking the Initial conditions of the Earth around the sun
P_X_E = 150e6;
P_Y_E = 100;
P_Z_E = 100;

V_X_E = 0.1;
V_Y_E = 30;
V_Z_E = 0.1;

Position_E = [P_X_E;P_Y_E;P_Z_E]; % form column vector for position
Velocity_E = [V_X_E;V_Y_E;V_Z_E]; % form column vector for velocity

%% Determining the orbital Mechanics

%Calculating angular momentum of the Earth in orbit km^2/s

he_vector = cross(Position_E,Velocity_E); %Angular momentum vector
h_e = norm(he_vector);     %Magnitude of the angular momentum

%We now find the eccentricity of the orbit

r_e = norm(Position_E)  %km
e_earth_vector = (cross(Velocity_E,he_vector)/u_sun) - (Position_E/r_e);
e_earth = 0.0167
% Semi-major axis of Earth's orbit

a_e = h_e^2/(u_sun*(1-e_earth^2))

%% Creating a state space model to simulate the Earth orbit numerically

%Initial system in ss form as there are 6 values that define the state of
%the system
sys_E = [Position_E(1);Position_E(2);Position_E(3);Velocity_E(1);Velocity_E(2);Velocity_E(3)];

no_of_orbits = 1;
% Function to handle the calculation of the time the body is in orbit
[T_E,Orbit_E] = Get_TOF(a_e,u_sun,e_earth,sun_Mass,Dist_to_Sun_self,no_of_orbits);

% Create a time interval
dt = 300; % time step
tspan_e = 0:dt:T_E;  %s

%% Propogation of Earth orbit

%Simulating using the ode45 solver 
options = odeset('RelTol',1e-8,'AbsTol',1e-9);
[t,S] = ode45(@(t,S) OrbitState(t,S,u_sun,sun_Radius),tspan_e,sys_E,options);
state = S;

%% Capturing the Orbit Trajectory

States_X_E = S(:,1); % all x
States_Y_E = S(:,2); % all y
States_Z_E = S(:,3); % all z
States_Xdot_E = S(:,4); % all xdot
States_Ydot_E = S(:,5); % all ydot
States_Zdot_E = S(:,6); % all zdot


% Running simulations for Orbit around the Earth

%% Taking the Initial conditions of the satellite around the Earth
switch Orbit_type
    case 'Moon_orbit'
        P_X = 400e3;
        P_Y = 100;
        P_Z = 100;

        V_X = 0.1;
        V_Y = 1.022;
        V_Z = 0.1;
    case 'Molniya_Orbit'
        P_X = 0;
        P_Y = -7000;
        P_Z = -14000;
% 
%         V_X = 2;
%         V_Y = 1;
%         V_Z = -0.5;
         V_X = 1.2035;
         V_Y = -6.175;
         V_Z = -0.30875;
        
    case 'Geo Stationary'
        P_X = 42241;
        P_Y = 0;
        P_Z = 0;

        V_X = 0;
        V_Y = 3.072;
        V_Z = 0;
    case 'Polar Orbit'
        P_X = 42241;
        P_Y = 0;
        P_Z = 0;

        V_X = 0;
        V_Y = 0;
        V_Z = 3.072;
     case 'Low Earth Orbit'
        P_X = 7000;
        P_Y = 0;
        P_Z = 0;

        V_X = 0;
        V_Y = 7.5/1.5;
        V_Z = 7.5/1.5;
        
    case 'High Earth Orbit'
        P_X = 3*7000;
        P_Y = 0;
        P_Z = 0;

        V_X = 0;
        V_Y = 7.5/1.5;
        V_Z = 7.5/1.5;
    otherwise
        %If orbit type is not recognized then use moon orbit type as the
        %default orbit parameters
        P_X = 150e6;
        P_Y = 100;
        P_Z = 100;

        V_X = 0.1;
        V_Y = 30;
        V_Z = 0.1;

end
Position = [P_X;P_Y;P_Z] % form column vector for position
Velocity = [V_X;V_Y;V_Z] % form column vector for velocity

%% Determining the orbital Mechanics

%Calculating angular momentum of the Earth in orbit km^2/s

h_vector = cross(Position,Velocity); %Angular momentum vector
h = norm(h_vector);     %Magnitude of the angular momentum

%We now find the eccentricity of the orbit

r = norm(Position)  %km
e_vector = (cross(Velocity,h_vector)/u_earth) - (Position/r);
e = norm(e_vector)
% Semi-major axis of Earth's orbit

a = h^2/(u_earth*(1-e^2))

%% Creating a state space model to simulate the Earth orbit numerically

%Initial system in ss form as there are 6 values that define the state of
%the system
sys_S = [Position(1);Position(2);Position(3);Velocity(1);Velocity(2);Velocity(3)];

no_of_orbits = 1;
% Function to handle the calculation of the time the body is in orbit
[T_S,Orbit] = Get_TOF(a,u_earth,e,earth_Mass,Dist_to_Sun,no_of_orbits);

% Create a time interval
dt = 300; % time step
tspan_S = 0:300:T_S;  %s

%% Propogation of Earth orbit

%Simulating using the ode45 solver 
options = odeset('RelTol',1e-8,'AbsTol',1e-9);
[t_S,S_S] = ode45(@(t_S,S_S) OrbitState(t_S,S_S,u_earth,earth_Radius),tspan_S,sys_S,options);
% [t,S] = ode45(@(t,S) OrbitState(t,S,u_sun,sun_Radius),tspan_e,sys_E,options);
state_S = S_S;

%% Capturing the Orbit Trajectory

States_X_S = state_S(:,1); % all x
States_Y_S = state_S(:,2); % all y
States_Z_S = state_S(:,3); % all z
States_Xdot_S = state_S(:,4); % all xdot
States_Ydot_S = state_S(:,5); % all ydot
States_Zdot_S = state_S(:,6); % all zdot


%% To delete all variables except for variables of interest run the code below

clearvars -except States_X_S States_Y_S States_Z_S States_X_E States_Y_E States_Z_E Orbit_type case
save('traj.mat')

%% 3D simulation
if Orbit_type ~= "Moon_orbit"
    Sun_Earth
elseif Orbit_type == "Moon_orbit"
    Sun_Earth_Moon
end
%% Functions for the simulation

%To get time of flight
function [T, Orbit] = Get_TOF(a,u,e,Body_Mass,Dist_to_Sun,no_of_orbits)

if e<1 && e>=0
    T = 2*pi*no_of_orbits*sqrt(a^3/u);
    if e==0
        Orbit = 'A circular orbit' ; 
    else
        Orbit = 'An elliptical Orbit';
    end
elseif e>1
    mass_of_sun = 1.989e30;
    R_soi = Dist_to_Sun*(Body_Mass/mass_of_sun)^(2/5);
    
    R_mag = R_soi;
    theta = acos(((a*(1-e^2)/R_mag)-1)/e);
    
    F = 2*atanh(sqrt((e-1)/(e+1))*tan(theta/2));
    T = sqrt((-a)^3/u)*(e*sinh(F) - F);
    Orbit = 'A Hyperbolic Orbit';
elseif e==1
    mass_of_sun = 1.989e30;
    R_soi = Dist_to_Sun*(Body_Mass/mass_of_sun)^(2/5);
    
    R_mag = R_soi;
    theta = acos(((a*(1-e^2)/R_mag)-1)/e);
    
    D = tan(theta/2);
    T = sqrt(2*(a*(1-e^2))^3/u)*(D + D^3/3);
    Orbit = 'A Parabolic Orbit';
else
    T = 0;
    Orbit = 'No Orbit';
    msgbox('Negative or complex eccentricity... Reenter data!','Data Error');
end
end

function [State] = OrbitState(t,S,u,Body_Radius)

r_new = norm(S(1:3));

A_d = [0.0;0.0;0.0];

State = [S(4);S(5);S(6);(-u/r_new^3)*S(1) + A_d(1);(-u/r_new^3)*S(2) + A_d(2); (-u/r_new^3)*S(3) + A_d(3)];
    
end





