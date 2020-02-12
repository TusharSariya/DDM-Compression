% map_in_space_new builds the map of the scattered power in space, given the transmitter and receiver positions...
% ...and velocities, and the MSS along the x and y direction, and returns the grid with the values of the relative...
% ...delays, the absolute doppler shifts and the doppler frequency relative to the specular point, and the grid which...
% ...accounts for the receiver antenna gain, the path losses and the bistatic radar cross section.

function [L_Kmeters,M_Kmeters,TotalXX,d_path,Dopp_off,c_dopp]=map_in_space_new(tx_pos,tx_vel,rx_pos,rx_vel,x)

%computes the DDMap in space

% Copywrite 2008, Scott Gleason & Maria Paola Clarizia
% license: GPL, see gpl.txt


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIAL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%GPS signal parameters
L1 = 1575.42e6;                     % GPS L1 carrier frequency [Hz]
lamda = 0.1904;                     % GPS L1 wavelength [m]

%geometrical parameters
WGS84_a = 6378137;		            % WGS84 earth radius [m]
r = 6356752.3142;	                % mean sea level [m]

%physical parameters
speedlight = 2.99792458*10^8;       % speed of light [m/s]
mu = 3.9860050e14;                  % used to compute example data

%temporal parameters
Integration_Time = 0.001;           % integration time [m]

sigmaL=x(1);
sigmaM=x(2);
%phi_wind=30;

% Calculate specular point location in WGS84 ECEF
[SP,Sxllh] = RT2S_Example1(rx_pos,tx_pos);
% rotate Rx and Tx to the model centered, surface based, reference frame
[Rx,Rxv,Tx,Txv,Sxangle_deg] = wgsxyz2model2(rx_pos,rx_vel,tx_pos,tx_vel,SP);  
Sxangle=d2r(Sxangle_deg);

Sx=[0 0 r];             % Specular Point position in model reference frame
Sxv=[0 0 0];            % Specular Point velocity

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRANSMITTER ANTENNA GAIN (CONSTANT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find the transmitter antenna gain (assumed constant)
% Tx2Sx_unit = (Tx - Sx)/norm(Tx - Sx);
% Tx_unit = Tx/norm(Tx);
% Tx_ang = abs(acos(dot(Tx2Sx_unit,Tx_unit)));                                     % Angle between the Tx-Sx line and the Transmitter position vector
% 
% TxG_table = [12.7 13 13.6 14.1 14.9 15.1 14.8 13.8 12.7 9.9 6.5 0.7];  
% X = [0 2 4 6 8 10 12 14 16 18 20 22];  
% [P] = polyfit(X,TxG_table,4);                                                    % 4th order curve fit
% 
% TxG = P(1)*(Tx_ang^4) + P(2)*(Tx_ang^3) + P(3)*(Tx_ang^2) + P(4)*Tx_ang + P(5);  % Antenna gain [dB], assumed constant over the patch 
% 
% TxG_temp = 10^(TxG/10)                                                          % convert to power in Watts                                                                                  

TxG_temp=1;  % or, keep it simple

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINITION AND GRIDDING OF THE OCEAN PATCH 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Width = 100000;                 % patch width [m]
%Height = 100000;                % patch height [m]
Width = 100000;                 % patch width [m]
Height = 100000;                % patch height [m]
Offset_width = 0;               % offset width [m] 
Offset_height = 0;              % offset height [m] 
resolution = 1000;	            % grid linear step [m]
Area_dS = resolution^2;         % grid element area [m^2]
magSx=norm(Sx);                 % Sx magnitude [m]

% use angles to account for the earth curvature

Width_angle = Width/magSx;                          % grid element angular width [rad]
Height_angle = Height/magSx;                        % grid element angular height [rad]
step_angle = resolution/magSx;	                    % grid angular step [rad]

Offset_width_angle = Offset_width/magSx;            % offset width angle [rad]
Offset_height_angle = Offset_height/magSx;          % offset height angle [rad]
                                                        
Instart = -Width_angle + Offset_width_angle;        % starting horizontal angle [rad]
Inend = Width_angle + Offset_width_angle;           % ending horizontal angle [rad]

Outstart = -Height_angle + Offset_height_angle;     % starting vertical angle [rad]
Outend = Height_angle + Offset_height_angle;        % ending vertical angle [rad]

In_meters = -Width + Offset_width;                  % patch width + offset [m]
Out_meters = -Height + Offset_height;               % patch height + offset [m] 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOMINAL DISTANCE TRANSMITTER-SPECULAR POINT-RECEIVER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nominal_distance = norm(Rx-Sx) + norm(Tx-Sx);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DOPPLER FREQUENCY AT THE SPECULAR POINT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xxx,c_dopp] = calcdopps(Rx,Rxv,Tx,Txv,Sx,Sxv);          % xxx = direct doppler frequency (not used)
                                                           % c_dopp = reflected doppler frequency
                                                            

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DOUBLE LOOP TO CALCULATE THE SCATTERING POWER FOR EACH GRID ELEMENT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

width_index = 1;     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% width loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = Instart:step_angle:Inend
    
%    i
    phi_rad=i;
    phi_deg(width_index) = r2d(i);              % express width angle in degrees
    
    L_Kmeters(width_index) = (magSx*phi_rad)/1000;   % linear horizontal distance between the specular point and the grid element [Km]

    length_index = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% length loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for j = Outstart:step_angle:Outend
        
            theta_rad=j;
            theta_deg(length_index) = r2d(j);       % express height angle in degrees
        
            
            M_Kmeters(length_index) = (magSx*theta_rad)/1000;   % linear vertical distance between the specular point and the grid element [Km]
        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % find the coordinates of the scattering point, using rotation matrices
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            T_II2x_A = [-cos(phi_rad) 0 sin(phi_rad) 	
		                0	1		0  
	            sin(phi_rad) 0 cos(phi_rad)];           % rotation with respect to the width angle
  
  
  
            T_II2x_B = [1 		0			0 	
	                0	cos(theta_rad) -sin(theta_rad) 	
	            0	sin(theta_rad) cos(theta_rad)];     % rotation with respect to the height angle
  

            PUT = Sx*T_II2x_A*T_II2x_B;                 % coordinates of the scattering point (grid element)
    
    
            RSx_unit = (Rx-PUT)./norm(Rx-PUT);           % transmitter-scattering point vector
            TSx_unit = (Tx-PUT)./norm(Tx-PUT);           % receiver-scattering point vector


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % calculation of Sigma0
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            q_vec=((2*pi*L1)/speedlight)*(TSx_unit+RSx_unit);       % scattering vector (bisector of the angle between the tx-scattering point line and the rx-scattering point line)
            [temp] = sigma0mp(q_vec,r2d(Sxangle),sigmaL,sigmaM);    % calculate the bistatic radar cross section (Scott's thesis,pp.18)
   
            sigma0(length_index,width_index) = temp;
            sigma0dS_temp = sigma0(length_index,width_index)*Area_dS;   %sigma0 over the grid element
    
    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Delays and Path Loss 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            R1(length_index,width_index) = norm(Tx-PUT);     % Tx-scattering point range
            R2(length_index,width_index) = norm(Rx-PUT);     % Rx-scattering point range

            Path_TxS(length_index,width_index) = 1/(R1(length_index,width_index)^2);  
            Path_RxS(length_index,width_index) = 1/(R2(length_index,width_index)^2); 

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Receiver Antenna Gain
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            RxG = 1;
            RxG_temp = 10^(RxG/10);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  Doppler Frequencies calculations
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            [xxx,D1] = calcdopps(Rx,Rxv,Tx,Txv,PUT,Sxv);    % calculate the direct and reflected doppler frequency for the scattering point considered      

            Dopp_off(length_index,width_index)=D1;

            Diff_Dopp = D1 - c_dopp;                          % difference between the calculated doppler frequency and the one corresponding to the specular point

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % combine terms
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            dS_temp = TxG_temp^2*Path_TxS(length_index,width_index)*sigma0dS_temp*Path_RxS(length_index,width_index)*RxG_temp;
            Total_power(length_index,width_index) = dS_temp;            % total scattered power from the grid element
            
            length_index = length_index + 1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end     % length Loop end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    width_index = width_index + 1;
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
end     % width loop end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   gain=calc_antenna_gain(Tx,Rx,Rxv);
   gain=1;
   TotalXX=Total_power.*(gain);
   
   [rr,ss]=size(R1);
   d_path=(R1+R2)-(nominal_distance*ones(rr,ss));                    % Relative ranges
   
   

%    
%   figure(1);clf
% zoom on
% hold on
% surf(L_Kmeters,M_Kmeters,(R1+R2)./speedlight)
% title('Delays')
% colorbar
% 
% 
% figure(2);clf
% zoom on
% hold on
% surf(L_Kmeters,M_Kmeters,Dopp_off)
% title('Doppler frequencies')
% colorbar
% 
% 
% figure(3);clf
% zoom on
% hold on
% surf(L_Kmeters,M_Kmeters,sigma0)
% title('sigma0')
% colorbar
% 
% 
% 
% figure(6);clf
% zoom on
% hold on
% surf(L_Kmeters,M_Kmeters,TotalXX)
% xlabel('Km')
% ylabel('Km')
% %axis('equal')
% title('Scattering Power with 3dB Antenna Ellipse')
% colorbar
% 
% 
% 
%     
%     figure(4);clf
%     zoom on
%     hold on
%     surf(L_Kmeters,M_Kmeters,RxG_grid)
%     xlabel('Km')
%     ylabel('Km')
%     title('Receiver Antenna Gain')
%     colorbar
%     
% 
% 
% 
%     figure(7);clf
%     zoom on
%     hold on
%     surf(L_Kmeters,M_Kmeters,Path_TxS.*Path_RxS)
%     xlabel('Km')
%     ylabel('Km')
%     axis('equal')
%     title('Path losses')
%     colorbar
    



 
