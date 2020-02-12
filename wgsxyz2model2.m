function [Rxmodel,Rxvmodel,Txmodel,Txvmodel,Iangle_deg] = wgsxyz2model2(Rx,Rxv,Tx,Txv,Sx);

%[pos,vel] = wgsxyz2model(Rx,Tx,Sx);
%
% Inputs:
%	Rx = Receiver position and velocity in ECEF
%	Tx = Transmitter position and velocity  in ECEF
%	Sx = Specular Point position and velocity in ECEF
%
% Outouts:	
%
%	Rxmodel = Rx position in model reflerence frame  
%	Txmodel = Txx position in model reflerence frame  
%	Rxvmodel = Rx velocity in model reflerence frame  
%	Txvmodel = Txx velocity in model reflerence frame  
%	Iangle = Incendent angle (degrees)
%
% Copywrite 2008, Scott Gleason & Maria Paola Clarizia
% license: GPL, see gpl.txt

% Convert Sx to lat lon
[Sxlat,Sxlon,Sxalt] = wgsxyz2lla(Sx(1:3));

% Convert Rx and Tx into local ENU frame centered at Sx 
Rxenu=wgsxyz2enu((Rx)', Sxlat, Sxlon, Sxalt);
Txenu=wgsxyz2enu((Tx)', Sxlat, Sxlon, Sxalt);

% Convert Rxv and Txv into local ENU frame centered at Sx 
Rxvenu=wgsxyz2enu((Sx+Rxv)', Sxlat, Sxlon, Sxalt);
Txvenu=wgsxyz2enu((Sx+Txv)', Sxlat, Sxlon, Sxalt);

% Position, LMN

% find incedence angle
SxRx_unit = (Rx-Sx)/norm(Rx-Sx);
Sx_unit = Sx/norm(Sx);
Iangle = abs(acos(dot(SxRx_unit,Sx_unit)));
Iangle_deg = r2d(Iangle);
magRx = norm(Rx);
magTx = norm(Tx);
magSx = norm(Sx);

% find earth centre angles

Rx_unit = Rx./magRx;
Tx_unit = Tx./magTx;
Sx_unit = Sx./magSx;

c2 = acos(dot(Rx_unit,Sx_unit));
c1 = acos(dot(Tx_unit,Sx_unit));

r2d(c1);
r2d(c2);

Rxmodel = [(magRx)*sin(c2) 0 (magRx)*cos(c2)];		% Receiver

Txmodel = [-1*(magTx)*sin(c1) 0 (magTx)*cos(c1)];		% GPS satellite

% Velocity, determine angle between vel in ECEF and RxTx in model
Rx2Tx_unit = (Rx_unit - Tx_unit)./norm(Rx_unit - Tx_unit);

Rx2Txenu=wgsxyz2enu((Sx+Rx2Tx_unit)', Sxlat, Sxlon, Sxalt);
Rx2Tx_en_angle = r2d(atan2(Rx2Txenu(1),Rx2Txenu(2)));

Rxvel_ecef_en_angle = r2d(atan2(Rxvenu(1),Rxvenu(2)));
Txvel_ecef_en_angle = r2d(atan2(Txvenu(1),Txvenu(2)));

Rxalpha = Rxvel_ecef_en_angle - Rx2Tx_en_angle;
Rxalpha_rad = d2r(Rxalpha);

Txalpha = Txvel_ecef_en_angle - Rx2Tx_en_angle;
Txalpha_rad = d2r(Txalpha);

Rxv_hor = sqrt((Rxvenu(1)^2 + Rxvenu(2)^2));
Rxvmodel = [Rxv_hor*cos(Rxalpha_rad) Rxv_hor*sin(Rxalpha_rad) Rxvenu(3)];

Txv_hor = sqrt((Txvenu(1)^2 + Txvenu(2)^2));
Txvmodel = [Txv_hor*cos(Txalpha_rad) Txv_hor*sin(Txalpha_rad) Txvenu(3)];

