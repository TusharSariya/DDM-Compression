
This is some example code, which we admit is far from perfect, for generating simulated delay doppler maps.

start by running generate_DDM.m in either Octave or MATLAB.
The initial conditions are set at the start of this file using the variables,

for the transmiter position and velocity
tx_pos = [-11178791.991294 -13160191.204988 20341528.127540];
tx_vel = [2523.258023  -361.592839  1163.748104];

for the receiver position and velocity
rx_pos = [-4069896.7033860330 -3583236.9637350840 4527639.2717581640];
rx_vel = [-4738.0742342063 -1796.2525689964 -5654.9952013657];

The initial bidirectional sea state mean square slopes
MSS_start = [0.01 0.02];

Initial positions were selected based on actual data collections.

For example,
Space receiver positions can be found in the Nav_RL.dat file, at given times (see SoftwareReceiver_Data/Dec7_2005_Land).
corresponding satellite positions can be found in IGS_RL_sv15.dat for example. 

simialar data is included for ice and ocean data in the Feb4_2005_Data_Ice and Nov16_2005_Data_Ocean, respectively.

Alternatively, you could generate data using the GNSS Simulator.

To experiment with differnet sea state mean square slopes, we recommend using the Elfouhaily wave model scripts, also included.
These scripts will predict bidirectonal mean square slopes under various conditions (applogies, but the script has not been tested in Octave yet)
See initial publication for details. (Journal of Geophysical Research, 102, No. C7, 15781-15796.)

Please make any enhances or fixes available to us, it would be much appreciated.