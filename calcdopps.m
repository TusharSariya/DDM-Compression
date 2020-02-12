function [directdopp,reflecteddopp] = calcdopps(Recef,Rvecef,Tecef,Tvecef,posspec,velspec)

% This function solves for the doppler frequencies of 
% the direct and reflected signals
%
% Copywrite 2008, Scott Gleason & Maria Paola Clarizia
% GPL, see gpl.txt

L1 = 1575.42e6;             %GPS carrier frequency
speedlight = 2.99792458e8;
mcps = 1.023e6;  

magTX = norm(Tecef);
magRX = norm(Recef);

% direct doppler
% radial velocity = (Rvecef - Tvecef)dot(unit_RxTx)

RxTxv = (Rvecef + Tvecef);
RxTx = (Tecef-Recef)./norm(Tecef-Recef);

directv = dot(RxTxv,RxTx);

directdopp = directv*L1/speedlight; 
%directdopp = directv*mcps/speedlight; 

% reflected doppler
% max relative radial velocity = (velRx - velspec)dot(unit_Rxspec) +
%		(velspec - velTx)dot(unit_specTx) */

Rxspecv = Rvecef - velspec;
temp1 = posspec - Recef;
unit_Rxspec = temp1/norm(temp1);

specTxv = Tvecef-velspec;
temp2 = posspec-Tecef;
unit_specTx = temp2/norm(temp2);

gamma = (pi - (acos(dot(-unit_Rxspec,-unit_specTx))) )/2;     %elevation angle

reflectedv1 = dot(specTxv,unit_specTx);
Sxdopp = reflectedv1*L1/speedlight; 

reflectedv2 = dot(Rxspecv,unit_Rxspec);

%reflectedv2 = dot(Rxspecv,unit_Rxspec) + dot(specTxv,unit_specTx);

Rxdopp = reflectedv2*(L1+Sxdopp)/speedlight; 

reflecteddopp = Sxdopp + Rxdopp;



