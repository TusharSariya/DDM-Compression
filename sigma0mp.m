function [forw_sigma0] = sigma0mp(q_vec,Sxangle,sigma2x,sigma2z)

% method=4;
% thetaX=11.284;
% thetaZ=85.835;
% gamma=41.70;
% wind_speed=10;
% Sxangle=r2d(0.09);
% sigma2x=0.0147;
% sigma2z=0.0094;
%function [forw_sigma0] = sigma0_bistatic_3D(method,thetaX,thetaZ,phiS,wind_speed,Sxangle,sigma2x,sigma2z)
%
% Inputs:
%	method: Wind speed is converted to sea slopes using one of the following
%		1) Cox and Monk
%		2) Elfouhaily
%		3) Ebuchi and Kizu
%		4) Use provided values
%
%	thetaX: angle between incedent ray projected into the xy plane and the zy plane
%		this is the component from the sigmaz mean wave slope
%	thetaY: angle between incedent ray projected into the xz plane and the xy plane
%		this is the component from the sigmax mean wave slope
%	gamma: angle between Incedent/reflected ray bisector with the local normal
%
%	wind_speed: m/s
%
%	Sxangle = the angle resulting in specular reflection (deg)
%+
%	sigma2x = mean square upwind slopes, calculated elsewhere, if method=4
%		Elfouhaily calculation takes a little time, better to do it 1 time outside
%
%	sigma2z = mean square cross wind slopes, calculated elsewhere, if method=4
%		Elfouhaily calculation takes a little time, better to do it 1 time outside
%
% Outputs: 
%	forw_sigma0 - bistatic forward scattering sigma0 in dB's (Barrick Model)
%
% Copywrite 2008, Scott Gleason & Maria Paola Clarizia
% license: GPL, see gpl.txt

% alter the angles for bistatic scatter
% 
meanx = 0;
meanz = 0;

%x = d2r(thetaX - Sxangle);
 x=-q_vec(1)/q_vec(3);           
 y=-q_vec(2)/q_vec(3);             
 V=[x y];                          


%M=[cos(phi) -sin(phi);sin(phi) cos(phi)]*[sigma2x 0;0 sigma2z]*[cos(phi) sin(phi);-sin(phi) cos(phi)]; 
%inv_M=inv(M);                                                                                          
%W=1/(2*pi*sqrt(det(M)))*exp((-1/2)*V*inv_M*V');  


%W = (1/sqrt(2*pi*sqrt(sigma2x)*sqrt(sigma2z)))*exp(temp1+temp2);
W = 1/(2*pi*sqrt(sigma2x*sigma2z))*exp(-((x^2/(2*sigma2x))+(y^2/(2*sigma2z))));    
%W = (1/(sqrt(2*pi)*sqrt(sigma2x)))*exp(temp1);

%temp3 = pi*((mag(q_vec))^4)/(q_vec(3)^4);
temp3 = ((norm(q_vec))^4)/(q_vec(3)^4);
%temp3 = 1;

% Fresnel coeff
% emmisivity, A. P. Stogryn

tempdeg = 10; % sae temp degC
er = (3.70886e4 - 8.2168e1*tempdeg)/(4.21854e2 + tempdeg);
e = 73;  % dielectric constant - sea water at ~19cm

n1 = 1;   % index of refraction - air
n2 = sqrt(e); % index of refraction - sea water

reflect_angle = d2r(Sxangle);  
refract_angle = asin((n1/n2)*sin(reflect_angle));

E_refl = -1*sin(reflect_angle - refract_angle)/ ...
	sin(reflect_angle + refract_angle);

R = E_refl^2;
%R = 1;

forw_sigma0 =  temp3*W*R;



