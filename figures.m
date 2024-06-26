function figures(r0, v0, rf, vf, rc, vc, model)
% This user defined function takes in 3 sets of state vectors in the ECI
% frame and outputs two graphical representations. The first state vector
% represents the initial orbit and position of the main satellite. The
% second state vector represents the propagated orbit and position of the
% main satellite. And the third state vector represents the orbit and
% position of another satellite to be compared to the main satellite.
% Therefore, 2 graphs are outputted.
%
% INPUTS GUIDE:
% format: figures(r0, v0, rf, vf, rc, vc, model)
% --> r: ECI position vector [km]
% --> v: ECI velocity vector [km/s]
%    --> subscript 0 denotes the initial states of the main satellite
%    --> subscript f denotes the propagated states of the main satellite
%    --> subscript c denotes the sates of the compare satellite
% --> model: char varaible specifying the Earth model used in calculating
% the initial and propagated state vectors of the main satellite; 'oblate'
% for Oblate-speroid model or 'sphere' for spherical model
%
% OUTPUTS GUIDE:
% --> This function will output two graphs. One graph containing the
% intitial and propagated positions and orbits of the main satellite, and
% one graph showing the propagated orbit and position of the main satellite
% compared to the orbit and position of the compare satellite.
% NOTE: If a spherical model was used for the main satellite, only one
% orbit will be displayed on the first graph as the orbit does not move in
% 3 dimensional space for that model.

% Define the Earth as a sphere for plotting
R = 6378.1366; % [km]
[xE, yE, zE] = sphere;
xE = xE * R; % [km]
yE = yE * R; % [km]
zE = zE * R; % [km]

% Convert state vectors into orbital elements using the rv2oe local
% function
oe0 = rv2oe(r0, v0);
oef = rv2oe(rf, vf);
oec = rv2oe(rc, vc);

% Define the x, y, and z arrays for the initial and propagated position
% vectors for subsequent plotting
x0 = linspace(0,r0(1), 361); % [km]
y0 = linspace(0,r0(2), 361); % [km]
z0 = linspace(0,r0(3), 361); % [km]
xf = linspace(0,rf(1), 361); % [km]
yf = linspace(0,rf(2), 361); % [km]
zf = linspace(0,rf(3), 361); % [km]
xc = linspace(0,rc(1), 361); % [km]
yc = linspace(0,rc(2), 361); % [km]
zc = linspace(0,rc(3), 361); % [km]

% Check to see if initial and final orbits are the same
if strcmp(model, 'sphere')
    % If the intitial and final elements are the same then the orbit has
    % not changed in three dimensional space; thus only plot one orbit with
    % two positions
    
    % Initialize array to hold position vectors
    r_vect_main = []; % [km]
    r_vect_comp = [];
    % Define a loop for one revolution to find the position at each value
    % of the true anomaly
    for f = 0:1:360 % [deg]

        % Call the oe2r function for each true anomaly value to find the
        % corresponding ECI positions
        oe0(6) = f; % [deg]
        r0_temp = oe2r(oe0); % [km]
        oec(6) = f; % [deg]
        rc_temp = oe2r(oec); % [km]

        
        % Save this iterations r value into the r_vect array for subsequent
        % plotting
        r_vect_main = [r_vect_main, r0_temp]; % [km]
        r_vect_comp = [r_vect_comp, rc_temp]; % [km]

    end
    
    % Plot the orbit and the intitial and final positions of main satellite
    figure
    plot3(r_vect_main(1, :), r_vect_main(2, :), r_vect_main(3, :), 'b', ...
        x0, y0, z0, 'k:', xf, yf, zf, 'm:', 'linewidth', 2.5), grid, ...
        xlabel('\bf I (km)'), ylabel('\bf J (km)'), zlabel('\bf K (km)'), ...
        title('ECI Spherical Model Orbit Plot')
    hold on
    surf(xE, yE, zE, 'EdgeColor', 'none')
    alpha 0.2
    legend('Orbit', 'Initial Position', 'Propagated Position', ...
        'location', 'best')
    hold off

    % Plot the propagated orbit and position of the main satellite with the
    % position and orbit of the compare satellite
    figure
    plot3(r_vect_main(1, :), r_vect_main(2, :), r_vect_main(3, :), 'b', ...
        r_vect_comp(1, :), r_vect_comp(2, :), r_vect_comp(3, :), ...
        xf, yf, zf, 'k:', xc, yc, zc, 'm:', 'linewidth', 2.5), grid, ...
        xlabel('\bf I (km)'), ylabel('\bf J (km)'), zlabel('\bf K (km)'), ...
        title('ECI Plot of Main and Compare Satellites at Propagated Epoch')
    hold on
    surf(xE, yE, zE, 'EdgeColor', 'none')
    alpha 0.2
    legend('Main Sat. Orbit', 'Compare Sat. Orbit', 'Main Sat. Pos.', ...
        'Compare Sat. Pos.','location', 'best')
    hold off


else
    % If initial and final elements are not the same then orbit has changed
    % in three dimensional space; thus plot two orbits and two positions
    
    % Initialize array to hold position vectors fo both the initial and
    % propagated states
    r0_vect_main = []; % [km]
    rf_vect_main = []; % [km]
    r_vect_comp = []; % [km]
    
    % Define a loop for one revolution to find the position at each value
    % of the true anomaly
    for f = 0:1:360 % [deg]

        % Call the oe2r function for each true anomaly for both the initial
        % and propagated states to find the ECI positions at each true
        % anomaly
        oe0(6) = f; % [deg]
        r0 = oe2r(oe0); % [km]
        oef(6) = f; % [deg]
        rf = oe2r(oef); % [km]
        oec(6) = f; % [deg]
        rc_temp = oe2r(oec); % [km]

        % Save this iterations value into the respective position arrays
        r0_vect_main = [r0_vect_main, r0]; % [km]
        rf_vect_main = [rf_vect_main, rf]; % [km]
        r_vect_comp = [r_vect_comp, rc_temp]; % [km]

    end

    % Plot the Initial and final orbits as well as the intitial and final
    % positions
    figure
    plot3(r0_vect_main(1, :), r0_vect_main(2, :), r0_vect_main(3, :), 'b', ...
        rf_vect_main(1, :), rf_vect_main(2, :), rf_vect_main(3, :), ...
        x0, y0, z0, 'k:', xf, yf, zf, 'm:', 'linewidth', 2.5), grid, ...
        xlabel('\bf I (km)'), ylabel('\bf J (km)'), zlabel('\bf K (km)'), ...
        title('ECI Oblate Spheroid Model Orbit Plot')
    hold on
    surf(xE, yE, zE, 'EdgeColor', 'none')
    alpha 0.2
    legend('Initial Orbit', 'Propagated Orbit', 'Initial Position', ...
        'Propagated Position', 'location', 'best')
    hold off

    % Plot the propagated orbit and position of the main satellite with the
    % position and orbit of the compare satellite
    figure
    plot3(r0_vect_main(1, :), r0_vect_main(2, :), r0_vect_main(3, :), 'b', ...
        r_vect_comp(1, :), r_vect_comp(2, :), r_vect_comp(3, :), ...
        xf, yf, zf, 'k:', xc, yc, zc, 'm:', 'linewidth', 2.5), grid, ...
        xlabel('\bf I (km)'), ylabel('\bf J (km)'), zlabel('\bf K (km)'), ...
        title('ECI Plot of Main and Compare Satellites at Propagated Epoch')
    hold on
    surf(xE, yE, zE, 'EdgeColor', 'none')
    alpha 0.2
    legend('Main Sat. Orbit', 'Compare Sat. Orbit', 'Main Sat. Pos.', ...
        'Compare Sat. Pos.','location', 'best')
    hold off

end
end



% Define the oe2rv local function to take in orbital elements and calculate 
% the corresponding state vector
function r = oe2r(oe)

% Use the standard Orbital Elements to ECI State Vector Procedure in order
% to find the propagated state vector
% Define useful parameters
p = oe(1)*(1 - oe(2)^2); % [km]
rmag = p/(1 + oe(2)*cosd(oe(6))); % [km]
r_PQW = [rmag*cosd(oe(6)); rmag*sind(oe(6)); 0]; % [km]


% Establish important variables for rotation matrix R
col = cosd(oe(5));
cou = cosd(oe(4));
ci = cosd(oe(3));
sol = sind(oe(5));
sou = sind(oe(4));
si = sind(oe(3));

% Define the rotation matrix R
R = [cou*col - sou*sol*ci, -cou*sol - sou*col*ci, sou*si; sou*col + cou*sol*ci, ...
    -sou*sol + cou*col*ci, -cou*si; sol*si, col*si, ci];

% Covert the PQW vector to ECI
r = R*r_PQW; % [km]

end



% Define the rv2oe local function to take in the state vector and calculate
% the corresponding orbital elements
function oe = rv2oe(r, v)
mu = 398600.4354; % [km^3/s^2]

% Find the Orbital Elements using standard ECI to KOE procedure
% Create empty vector for keplerian elements
oe = zeros(6, 1);

% Coordinate frame unit vectors
I = [1 0 0]; J = [0 1 0]; K = [0 0 1];

rhat = r/norm(r); % position unit vector [km]
h = cross(r,v); % angular momentum {h = r x v}
hhat = h/norm(h); % normalized angular momentum
nhat = cross(K,h)/norm(cross(K,h)); % normalized ascending node vector

% Eccentricity (e)
e = (1/mu)*cross(v,h) - rhat;
oe(2) = norm(e);
energy = (1/2)*dot(v,v) - mu/norm(r); % [km^2/s^2]

% Semi-major axis (a)
oe(1) = -mu/(2*energy); % [km]

% Inclination (i) of orbit
oe(3) = acos(dot(K,hhat)); % [rad]

% Right ascension of the ascending node (Omega)
oe(4) = mod(atan2(dot(J,nhat),dot(I,nhat)), 2*pi); % [rad]

% Argument of periapsis (w)
oe(5) = mod(atan2(dot(hhat,cross(nhat,e)),dot(nhat,e)), 2*pi); % [rad]

% True anomaly (f) at epoch [rad]
oe(6) = mod(atan2(dot(hhat,cross(e,r)),dot(e,r)), 2*pi); % [rad]

% Convert all angles to degrees for final output
oe(3:end) = rad2deg(oe(3:end)); % [deg]


end
