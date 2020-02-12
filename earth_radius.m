function r = earth_radius(ECEF_pos)

%function r = earth_radius(ECEF_pos)
%
%	ECEF_pos = position vector in ECEF 
%	r = radius of earth on that line, in meters
%
% Copywrite 2008, Scott Gleason
% GPL, see gpl.txt

WGS84_a = 6378137;
WGS84_f = 1/298.257223563;
WGS84_b = WGS84_a*(1-WGS84_f);
WGS84_e = sqrt(2*WGS84_f - WGS84_f^2);

% calculate latitude angle
magpos = norm(ECEF_pos);
z = ECEF_pos(3);  % z axis
theta = asin(z/magpos);

% calculate ellipse radius
temp1 = 1 - WGS84_e^2;
temp2 = 1 - ((WGS84_e^2)*(cos(theta)^2));
r = WGS84_a*sqrt(temp1/temp2);


