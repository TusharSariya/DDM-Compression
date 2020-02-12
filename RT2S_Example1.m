function [S,S_llh] = RT2S_Example1(R,T)
%
% function [S] = RT2S_Example1(R,T)
%
% Input: 
%	R = Receiver Position	ECEF
%	T = GPS Satellite Position ECEF
% Output:
%	S = specular point position ECEF
%	S_llh = specular point position lat/lon/hgt
%
%
% Copywrite 2008, Scott Gleason
% license: GPL, see gpl.txt

close all

% Initial specular point guess, on the surface directly below R
r = earth_radius(R);
Rmag = norm(R);
Stemp = R*(r/Rmag);
S = Stemp;

% misc 
correction = 10000;
tol = 0.001; 		% convergence tolerance, meters
iterations = 1;		
rad2deg = 180/pi;
K = 10000;		% correction gain

while correction > tol

	% Take derivatives
	S2R_unit = (R - S)./norm(R-S);
	S2T_unit = (T - S)./norm(T-S);
	% calculate correction direction
	dXYZ = S2R_unit + S2T_unit;

	% apply raw correction
	Stemp = S + K*dXYZ;
	
	% contrain to Earth surface
	r = earth_radius(Stemp);
	Stemp = (Stemp./norm(Stemp))*r;

	% watch real correction magnitudes, should continually decrease
	correction_temp = abs(norm(S-Stemp));
	correction(iterations) = correction_temp;

	% update estimate of specular point
	S = Stemp;

	% adjust gain based on correction, i.e. if we are getting closer, lower correction gain
	if(correction_temp > 10)
		K = 10000;
	else
		K = 1000;
	end

	iterations = iterations + 1;
	if(iterations > 10000)
		break;
	end

end

% convert final value to lat lon hgt
S_llh = [0 0 0]';
[S_llh(1),S_llh(2),S_llh(3)] = wgsxyz2lla(S);


% extra

% sanity check, test Snell's Law
RS_unit = (R-S)./norm(R-S);
TS_unit = (T-S)./norm(T-S); 
S_unit = S./norm(S);
theta1 = rad2deg*(acos(dot(RS_unit,S_unit)));
theta2 = rad2deg*(acos(dot(TS_unit,S_unit)));
theta_diff_degrees = abs(theta1 - theta2);

%figure(200)
%plot(correction)
%xlabel('iteration')
%ylabel('correction, meters')
%title('Correction Convergence')
%zoom on
