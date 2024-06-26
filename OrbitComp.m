
function [r0, v0, oe0, rf, vf, oef] = OrbitComp(lat, lst, alt, ra, dec, JD, JD_prop, model)
% The following user-defined function employs Gauss' method of preliminary
% orbit determination to determine the orbital elements and state vectors
% at an initial observation and at a propagated position. To do this, the
% function is given the observation site lattitude, 3 local siderial times
% corresponding to the 3 observation times, the site altitude, and 3 sets
% each of the right ascension, declination, and julian dates corresponding
% to each observation. In addition, the julian date corresponding to the
% desired propagated position is also inputted. Lastly, the type of model
% to be used is inputted as either 'oblate' for an oblate earth model or as
% 'sphere' for a sphereical earth model to be used in subsequent
% calculations. This is explained below in the function outline.
%
%
% NOTE: This function takes into account Earth's oblateness when
% calculating the site vectors as well as its cause on orbital precession
% (secular changes of the longitude of ascending node, argument of
% periapsis, and the mean/true anomaly) when model is set to 'oblate'
% Otherwise, model is set to 'sphere' and the site vectors are calculated
% using a spherical Earth model and secular changes for the longitude of
% ascending node, argument of periapsis, and mean anomaly are not taken
% into account.
%
% FUNCTION OUTLINE:
% 1) Use Gauss' angles only method to determine the 3 position vectors
% corresponding to each observation
% --> If model set to 'oblate' the site vectors will take into account
% Earth's oblateness, otherwise if model is set to 'sphere' site vectors
% will assume a sphereical model of Earth
% 2) Use either Gibb's or Herrick-Gibb's method to determine the state
% vector at the second observation
% 3) Convert the observation 2 state vector into the initial orbital
% elements
% 4) Determine the changes in the longitude of ascending node and argument
% of periapsis due to earths oblateness and use them to determine their
% propagated counterparts if model variable is set to 'oblate'
% --> If model set to 'sphere' secular changes for the above orbital
% elements will not be calculated
% 5) Use a Kepler's problem algorithm to find the propagated true anomaly
% according to both the mean motion and the secular mean anomaly change due
% to Earth's oblateness
% 6) Convert the propagated orbital elements into the corresponding
% propagated state vector
%
%
% INPUTS GUIDE: 
% --> lat: site latitude [deg]
% --> lst: 3x1 array containing the site local siderial times [deg]
% --> alt: site altitude [m]
% --> ra: 3x1 array containing the three right ascensions [deg]
% --> dec: 3x1 array containing the three declinations [deg]
% --> JD: 3x1 array containing the three julian dates
% --> JD_prop: propagated julian date
% --> model: string containing either 'oblate' or 'sphere'
%
% OUTPUTS GUIDE:
% --> r0: 3x1 position vector corresponding to LOS 2 [km]
% --> v0: 3x1 velocity vector corresponding to LOS 2 [km/s]
% --> oe0: 6x1 array containing LOS 2 orbital elements [UNITS VARY]
%     --> oe0(1): semimajor axis, a [km]
%     --> oe0(2): eccentricity, e [N/A]
%     --> oe0(3): inclination, i [deg]
%     --> oe0(4): longitude of the ascending node, cap_omega [deg]
%     --> oe0(5): argument of periapsis, low_omega [deg]
%     --> oe0(6): true anomaly, f [deg]
%
% --> rf: 3x1 position vector corresponding to progated position [km]
% --> vf: 3x1 velocity vector corresponding to propagated position [km/s]
% --> oef: 6x1 array containing propagated orbital elements [UNITS VARY]
%     --> oef(1): semimajor axis, a [km]
%     --> oef(2): eccentricity, e [N/A]
%     --> oef(3): inclination, i [deg]
%     --> oef(4): longitude of the ascending node, cap_omega [deg]
%     --> oef(5): argument of periapsis, low_omega [deg]
%     --> oef(6): true anomaly, f [deg]




% Define planetary constants
mu = 398600.4354; % [km^3/s^2]
R = 6378.1366; % [km]
J2 = 1.0826359e-3;

% Determine flattening factor depending on which model is desired
f = 0;
if strcmp(model, 'oblate')
    f = 0.003353;
else
    f = 0;
end

% Convert the altitude of the site in kilometers
hsite = alt/1e3; % [km]




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          GAUSS' METHOD                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Given the Julian dates, determine the change in times between the first
% and second observations, and the second and third observations
% Convert julian dates to seconds
ut = JD.*86400; % [sec]
% Determine change in times
T1 = ut(1, 1) - ut(2, 1); % [sec]
T3 = ut(3, 1) - ut(2, 1); % [sec]


% Determine the change in time helpful parameters
a1 = T3/(T3-T1); a3 = -T1/(T3-T1); % [sec]
a1u = (T3*((T3 - T1)^2 - T3^2))/(6*(T3 - T1)); % [sec^2]
a3u = (-T1*((T3 - T1)^2 - T1^2))/(6*(T3 - T1)); % [sec^2]


% Systematically determine the LOS vectors using the right ascension and
% declination values and input into matrix L
L = nan(3, 1);
for k = 1:3
    Lk = [cosd(dec(k, 1))*cosd(ra(k, 1));
        cosd(dec(k, 1))*sind(ra(k, 1));
        sind(dec(k, 1))];
    L = cat(2, L, Lk);
end
L = L(:,2:end);
% Determine the inverse LOS matrix
invL = (L)^-1;

% Determine the site position matrix according to the ellipsoidal earth
% model
r_site = nan(3, 1);
coeff1 = (R/sqrt(1 - (2*f - f^2)*sind(lat)^2) + hsite)*cosd(lat); % [km]
coeff2 = (R*(1 - f)^2)/(sqrt(1 - (2*f - f^2)*sind(lat)^2)) + hsite; % [km]
for k = 1:3
    r_sitek = [coeff1*cosd(lst(k)); coeff1*sind(lst(k)); coeff2*sind(lat)]; % [km]
    r_site = cat(2, r_site, r_sitek);
end
r_site = r_site(:, 2:end);

% Determine intermediate matrix M
M = invL*r_site;

% Define parameters for use in polynomial
d1 = M(2, 1)*a1 - M(2, 2) + M(2, 3)*a3; % [km]
d2 = M(2, 1)*a1u + M(2, 3)*a3u; % [km sec^2]
C = dot(L(:, 2), r_site(:, 2)); % [km]


% Define the 8th degree polynomial
p = [1, 0, -(d1^2 + 2*C*d1 + norm(r_site(:, 2))^2), 0, 0, ...
    -2*mu*(C*d2 + d1*d2), 0, 0, -mu^2*d2^2];
% Determine its roots and only save the real positive root
rts = roots(p); r2 = rts(real(rts) > 0 & imag(rts) == 0);


% calculate the value u and the c values
u = mu/(r2^3); 
c = [a1 + a1u*u; -1; a3 + a3u*u];


% Calculate the slant ranges
rhoc = M*-c;
rho0 = rhoc./c; % [km]

% Define the refinement iteration and count variable
% NOTE: loop terminates after 100 iterations for performance
k = 0;
while k < 100
    % Update the count index
    k = k + 1 ;

    % Determine the 3 position vectors in the ECI frame
    r1 = rho0(1, 1)*L(:, 1) + r_site(:, 1); % [km]
    r2 = rho0(2, 1)*L(:, 2) + r_site(:, 2); % [km]
    r3 = rho0(3, 1)*L(:, 3) + r_site(:, 3); % [km]
    
    % Call the gibbs function to find the state vector
    [r0, v0] = gibbs(r1, r2, r3, ut);
    
    % Call the rv2koe function to determine the orbital elements
    oe0 = rv2koe(r0, v0);
    
    % Define the parameter p
    p = oe0(1)*(1 - oe0(2)^2); % [km]

    % Define the lagrange coefficients for the angular difference between
    % r1 and r2 and r3 and r2
    ang1 = acos(dot(r1, r2)/(norm(r1)*norm(r2))); % [rad]
    ang3 = acos(dot(r3, r2)/(norm(r3)*norm(r2))); % [rad]
    f1 = 1 - (norm(r1)/p)*(1 - cos(ang1));
    f3 = 1 - (norm(r3)/p)*(1 - cos(ang3));
    g1 = norm(r1)*norm(r2)*sin(-ang1)/sqrt(mu*p);
    g3 = norm(r3)*norm(r2)*sin(ang3)/sqrt(mu*p);

    % Define new c values based off coefficients
    c(1) = g3/(f1*g3 - f3*g1);
    c(2) = -1;
    c(3) = -g1/(f1*g3 - f3*g1);

    % Calculate the slant ranges
    rhoc = M*-c;
    rho_test = rhoc./c; % [km]
 
    % Check if new slant ranges are equal to old slant ranges
    if isequal(round(rho_test, 10), round(rho0, 10))
            k = 101;
    end
    
    % Assign new slant ranges to initial value for another iteration if
    % needed
    rho0 = rho_test;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      PROPAGATED ORBITAL ELEMENTS                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the propagated orbital elements according to a non spherivcal
% earth model
% Determine the TOF
tof = JD_prop*86400 - ut(2, 1); % [sec]

% Initiate oef variable
oef = zeros(6, 1);

% Determine the final orbital elements
% --> for f = 0, spherical model is used and only the true anomaly changes
if f == 0
    oef(1:5, 1) = oe0(1:5, 1);
% --> for f = 0.003353, oblate model is used and there are secular changes
% in the ongitude of ascending node, argument of periapsis, and the mean
% anomaly
else
    oef(1:3, 1) = oe0(1:3, 1);
    
    % Determine the secular rates of the longitude of ascending node and
    % argument of periapsis according to an oblate Earth pertubation
    n = sqrt(mu/oef(1)^3); % angular rate [rad/s]
    wdot = 3/4*J2*n*(R/oef(1))^2*((5*cosd(oef(3))^2 - 1)/((1 - oef(2)^2)^2)); % [rad/s]
    odot = -3/2*J2*n*(R/oef(1))^2*(cosd(oef(3))/(1 - oef(2)^2)^2); % [rad/s]
    
    % Use the secular changes to determine the argument of periapsis and
    % longitude of scending node at the propogated julian date
    of = deg2rad(oe0(4)) + odot*(tof); % [rad]
    wf = deg2rad(oe0(5)) + wdot*(tof); % [rad]
    
    % Ensure values are between 0 and 360 degrees and define final output
    % variables
    oef(4) = mod(rad2deg(of), 360); % [deg]
    oef(5) = mod(rad2deg(wf), 360); % [deg]
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          KEPLERS PROBLEM                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Use the basis of keplers problem to propagate the orbit to the position
% marked by the JD input and find the new true anomaly.

% Define useful parameters
n = sqrt(mu/oe0(1)^3); % angular rate [rad/s]

% Find the initial and propogated eccentric and subsequent mean anomalies
E0 = 2*atan(sqrt((1 - oe0(2))/(1 + oe0(2)))*tand(oe0(6)/2)); % initial 
% eccentric anomaly [rad]
M0 = E0 - oe0(2)*sin(E0); % initial mean anomaly [rad]

% Define the secular initial mean anomaly change
M0dot = 3/4*J2*n*(R/oef(1))^2*((sqrt(1 - oef(2)^2)*(3*cosd(oef(3))^2 - 1))/ ...
    ((1 - oef(2)^2)^2)); % [rad/s]

% For f = 0, sphere model is used and the mean anomaly only changes based
% off the mean motion rate (n)
Mf = 0; % [rad]
if f == 0
    Mf = M0 + tof*n; % [rad]

% For f = 0.003353, oblate model is used and the mean anomaly changes based
% off of the secular change of the mean anomaly (M0dot) and the mean motion
% rate (n)
else
    Mf = M0 + tof*(n + M0dot); % [rad]
end

% Define the iteration index, starting guess, and the tolerance for
% convergence
k = 1;
Ef_0 = Mf + oe0(2)*sin(Mf) + oe0(2)^2/2*sin(2*Mf);
epsilon = 1e-8;

% Define the loop to iteratively solve for the root
Ef_k = Ef_0;
while abs(Ef_k - oe0(2)*sin(Ef_k) - Mf) > epsilon
    Ef_k = Ef_k - (Ef_k - oe0(2)*sin(Ef_k) - Mf)/(1 - oe0(2)*cos(Ef_k));
    k = k + 1;
end
Ef = Ef_k; % estimated eccentric anomaly for the propogated position [rad]

% Define the true anomaly of the propogated position
ff = mod(2*atand(sqrt((1 + oe0(2))/(1 - oe0(2)))*tan(Ef/2)), 360); % propogated 
% true anomaly [deg]


oef(6) = ff; % [deg]

% Call the koe2rv function to transform the propagated orbital elements
% into the correspinding state vector
[rf, vf] = koe2rv(oef); % [rf = km] ; [vf = km/s]


end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       PROBAGATED STATE VECTOR                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [r, v] = koe2rv(oe)
mu = 398600.4354; % [km^3/s^2]
% Use the standard Orbital Elements to ECI State Vector Procedure in order
% to find the propagated state vector
% Define useful parameters
p = oe(1)*(1 - oe(2)^2); % km
h = sqrt(p*mu); % km^2/s
rmag = p/(1 + oe(2)*cosd(oe(6))); % km
r_PQW = [rmag*cosd(oe(6)); rmag*sind(oe(6)); 0]; % km
v_PQW = [-mu/h*sind(oe(6)); mu/h*(oe(2) + cosd(oe(6))); 0]; % km/s

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

% Covert the PQW vectors to ECI
r = R*r_PQW; % km
v = R*v_PQW; % km/s

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             GIBBS METHOD                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [r0, v0] = gibbs(r1, r2, r3, ut)
mu = 398600.4354; % [km^3/s^2]
% Employ the GIBBS method to determine the velocity at r2

% Determine if the position vectors are coplanar
epsilon = (dot(r1, cross(r2, r3)))./(norm(r1)*norm(cross(r2, r3)));

% Find the angles between the position vectors
theta_12 = acosd((dot(r1, r2))./(norm(r1)*norm(r2))); % deg
theta_23 = acosd((dot(r2, r3))./(norm(r2)*norm(r3))); % deg

% Determine if the Gibbs Method Requirements for coplanarity and angle
% differences are met
% If they are, proceed into the gibbs method
v2 = 0;
if (abs(epsilon) < 0.0349) && (theta_12 > 1) && (theta_23 > 1)
    % Define the auxiliary vectors D, N, &  S
    D = cross(r2, r3) + cross(r3, r1) + cross(r1, r2);
    N = norm(r1).*cross(r2, r3) + norm(r2).*cross(r3, r1) + ...
        norm(r3).*cross(r1, r2);
    S = (norm(r2) - norm(r3)).*r1 + (norm(r3) - norm(r1)).*r2 + ...
        (norm(r1) - norm(r2)).*r3;

    % Find the ECI velocity vector at r2
    v2 = sqrt(mu/(norm(N)*norm(D))).*(1/norm(r2).*cross(D, r2) + S); % km
else
% Employ the Herrick-Gibbs Method if the Gibbs Method cannot be used
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       HERRICK-GIBBS METHOD                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Determine the time changes between positions
    t31 = ut(3, 1) - ut(1, 1); % [sec]
    t32 = ut(3, 1) - ut(2, 1); % [sec]
    t21 = ut(2, 1) - ut(1, 1); % [sec]
    v2 = -t32*(1/(t21*t31) + mu/(12*norm(r1)^3)).*r1 + ...
        (t32 - t21)*(1/(t21*t32) + mu/(12*norm(r2)^3)).*r2 + ...
        t21*(1/(t32*t31) + mu/(12*norm(r3)^3)).*r3; % [km/s]
end

r0 = r2; % [km]
v0 = v2; % [km/s]
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          ORBITAL ELEMENTS                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function oe0 = rv2koe(r0, v0)
mu = 398600.4354; % [km^3/s^2]
% Find the Orbital Elements using standard ECI to KOE procedure
% Create empty vector for keplerian elements
oe0 = zeros(6, 1);

% Coordinate frame unit vectors
I = [1 0 0]; J = [0 1 0]; K = [0 0 1];

rhat = r0/norm(r0); % position unit vector [km]
h = cross(r0,v0); % angular momentum {h = r x v}
hhat = h/norm(h); % normalized angular momentum
nhat = cross(K,h)/norm(cross(K,h)); % normalized ascending node vector

% Eccentricity (e)
e = (1/mu)*cross(v0,h) - rhat;
oe0(2) = norm(e);
energy = (1/2)*dot(v0,v0) - mu/norm(r0); % [km^2/s^2]

% Semi-major axis (a)
oe0(1) = -mu/(2*energy); % [km]

% Inclination (i) of orbit
oe0(3) = acos(dot(K,hhat)); % [rad]

% Right ascension of the ascending node (Omega)
oe0(4) = mod(atan2(dot(J,nhat),dot(I,nhat)), 2*pi); % [rad]

% Argument of periapsis (w)
oe0(5) = mod(atan2(dot(hhat,cross(nhat,e)),dot(nhat,e)), 2*pi); % [rad]

% True anomaly (f) at epoch [rad]
oe0(6) = mod(atan2(dot(hhat,cross(e,r0)),dot(e,r0)), 2*pi); % [rad]

% Convert all angles to degrees for final output
oe0(3:end) = rad2deg(oe0(3:end)); % [deg]

end