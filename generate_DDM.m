%
% Creates a DDM model, using the transmitter and receiver positions defined below.
%
% Copywrite 2008, Scott Gleason & Maria Paola Clarizia
% license: GPL, see gpl.txt

tx_pos = [-11178791.991294 -13160191.204988 20341528.127540];
tx_vel = [2523.258023  -361.592839  1163.748104];

rx_pos = [-4069896.7033860330 -3583236.9637350840 4527639.2717581640];
rx_vel = [-4738.0742342063 -1796.2525689964 -5654.9952013657];
%inc value for roughness of surface
MSS_start = [0.001 0.002];

%compute the map in spatial coordinates
[L,M,Ocean,d_path,dopptemp,centerdopp]=map_in_space_new(tx_pos,tx_vel,rx_pos,rx_vel,MSS_start);

%convert the map in delay-doppler coordinates
[Xt,Yt,Zt]=map_in_DD(Ocean,d_path,dopptemp,centerdopp);                

%compressed version
Xt_ = Xt(1:20:end);
Yt_ = Yt(5:5:end);
Zt_ = Zt(5:5:end,1:20:end);


% Plot 3D Delay Doppler Map
figure;
hold on
%subplot(2,2,1)
surf(Xt,Yt,Zt);
view(90,90);
axis([-inf inf 4000 inf])
title('Delay Doppler Map 1')
xlabel('Delay, (GPS L1 chips)*10')
ylabel('Doppler, Hz')

figure;
hold on
%subplot(2,2,1)
surf(Xt_,Yt_,Zt_);
view(90,90);
title('Delay Doppler Map 2')
xlabel('Delay, (GPS L1 chips)*10')
ylabel('Doppler, Hz')

