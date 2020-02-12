% map_in_DD creates the delay-doppler map, using information from map_in_space_new

function [X,Y,Znorm]=map_in_DD(Ocean,d_path,dopptemp,centerdopp)

% Copywrite 2008, Scott Gleason & Maria Paola Clarizia
% license: GPL, see gpl.txt

%close all
%clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIAL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LOcenterfreq = 1.405396825e6;  % approximate signal center frequency after the receiver processing
speedlight = 2.99792458*10^8;   % speed of light [m]
sealevel=6356752.3142;          % mean sea level [m]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%VARIABLES FROM FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

freq_step=5000;
dopp_bin_size=100;

max_samples = 3000;                         % maximum number of samples
sampleoffset = 50;                          % number of sample offset (at the beginning)

maxsampleindex = 1;                         % index of the maximum chip sample
samples_per_chip=5.58;

num_dopp_bins = 2*freq_step/dopp_bin_size;  % number of Doppler bins  
maxdoppindex = 1;                           % index of the maximum doppler shift found from data
mindoppindex = num_dopp_bins;               % index of the minimum doppler shift found from data

DD1_save = zeros(num_dopp_bins,max_samples+sampleoffset);   % data matrix (used to plot the DDmap)
DD1_count = zeros(num_dopp_bins,max_samples+sampleoffset);  % counter matrix (used to count how many delays in chips are binned per each doppler bin)

tempsize=size(Ocean);
XXstart = 1;
XXend = tempsize(1);
YYstart = 1;
YYend = tempsize(2);

chipindex = 1;

for XX = XXstart:XXend              % loop through Ocean Grid, X


    for YY = YYstart:YYend          % loop through Ocean Grid, Y

        
        if dopptemp(XX,YY)>centerdopp-freq_step & dopptemp(XX,YY)<centerdopp+freq_step 	    %choose only the frequency shifts contained in the chosen range 
            
	    value_1 = Ocean(XX,YY);

            %%%%DELAY BINNING%%%%
        
            delta_path = d_path(XX,YY);     % path range [m] (d_path contains the relative delays)
            del=delta_path/speedlight;      % delay in seconds
            del=del/0.97e-6;                % delay in chips (1 chip=1 microsec.)
            XYchip=del*samples_per_chip;        %total samples
     
 	    XYchip_index = floor(XYchip) + sampleoffset;    % find the column index in DD1_save corresponding to this delay

	    if(XYchip_index > maxsampleindex)
		    maxsampleindex = XYchip_index;              % update the index of the maximum chip sample
	    end

	    %%%%DOPPLER BINNING%%%%
      
	    dopp_offset = centerdopp -dopptemp(XX,YY);              % doppler shift (relative to the specular point)
	    XYdopp_bin = dopp_offset/dopp_bin_size;            
	    XYdopp_index = floor(XYdopp_bin) + num_dopp_bins/2+1;   % find the row index in DD1_save corresponding to this doppler shift

	    if(XYdopp_index > maxdoppindex)
		    maxdoppindex = XYdopp_index;                        % update the index of the maximum doppler shift
	    end
	    if(XYdopp_index < mindoppindex)
		    mindoppindex = XYdopp_index;                        % update the index of the minimum doppler shift
	    end
        
	    DD1_save(XYdopp_index,6+XYchip_index) = DD1_save(XYdopp_index,6+XYchip_index) + value_1;    % save the value of the scattered power in the correct delay index (column) and doppler bin (row)
	    DD1_count(XYdopp_index,6+XYchip_index) = DD1_count(XYdopp_index,6+XYchip_index) + 1;        % update the counter matrix
        
            DD1_save(XYdopp_index,4) = (floor(XYdopp_bin)*dopp_bin_size);     % save the value of the doppler bin TOLTO centerdopp-

    	    chipindex = chipindex + 1;

        end     % end if
        
	end     % end vertical scanning

end             % end horizontal scanning 

tempsize = length([mindoppindex:maxdoppindex])     % number of doppler bins found

%%%%%%%%%%NEW DATA MATRIX%%%%%%%%%%%%%

%DD2_save(:,1:5) = DD1_save(mindoppindex:maxdoppindex,1:5);                                   
DD2_save(:,6:5+maxsampleindex) = DD1_save(mindoppindex:maxdoppindex,6:5+maxsampleindex);    % hold the doppler bins for which there is scattered power

% Fill in the other columns
 
DD2_save(:,1) = 1;                      % sv
DD2_save(:,2) = 1;                      % start samples 
DD2_save(:,3) = maxsampleindex;         % end samples 
DD2_save(:,4) = DD1_save(:,4);
DD2_save(:,5) = 1;                      % ms 

%tempsize = size(DD2_save);
DD3_save=DD2_save;
tempsize=size(DD3_save);

Ti=0.001;           %apply to 10 samples before and after (1000/100)
           
%%%%%%%%%%%%%%%%%APPLY SINC FUNCTION TO EACH ELEMENT (1st version)%%%%%%%%%%%%%%%%%
         
        for l=6:tempsize(2) %do it for each delay
            
            for k=1:tempsize(1)
                
                bin_start=k-10;
                if bin_start<1
                    bin_start=1;
                end
                bin_end=k+10;
                if bin_end>tempsize(1)
                    bin_end=tempsize(1);
                end
                    
                bin_freq=DD2_save(k,4);
                   
                bin_start_freq=DD2_save(bin_start,4);
                bin_end_freq=DD2_save(bin_end,4);
                    
                x=[bin_start_freq:100:bin_end_freq];
                y=sinc((x-bin_freq).*Ti).^2;
                   
                value=y*DD2_save(bin_start:bin_end,l);
			
                DD3_save(k,l)=value;



            end
         end

                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
   DD4_save=DD3_save;
   tempsize=size(DD4_save);
   
%%%%%%%%%%%%%%%%%%%%%%%%APPLY CORRELATION FUNCTION(2nd version)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                
  for k=1:tempsize(1) %do it for each frequency bin
                 
                 for l=6:tempsize(2)
                 
                        
                         chip_start=l-5;
                         if chip_start<6
                             chip_start=6;
                         end
                         chip_end=l+5;
                         if chip_end>tempsize(2)
                             chip_end=tempsize(2);
                         end
                        
                         x=chip_start:chip_end;
                         vec_corr=lambda((x-l)./5.58).^2;
                         value=vec_corr*DD3_save(k,chip_start:chip_end)';
                         
                         DD4_save(k,l)=value;
                     end
                 end
                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sample_low=DD4_save(1,2);     
sample_high=DD4_save(1,3);     

low=6;                           
high=5+maxsampleindex;

X=[sample_low-sampleoffset:sample_high-sampleoffset].*1.79;                     % delay axis      
Y=DD4_save(:,4)'+abs(centerdopp);                                               % doppler axis
Z=DD4_save(:,low:high);                                                         % corresponding scattered power values

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%REMOVE NOISE%%%%%%

tempsize=size(Z);
[I,J]=find(Z);
noise_start = 1;
noise_end=J(1);             % find the first matrix column where there's at least one non-zero value

 meantemp = mean(Z(:,noise_start:noise_end),2);
% meantempmatr=meantemp*ones(1,tempsize(2));
% Znorm= Z- meantempmatr;                         % remove noise
%Z(:,noise_start:noise_end)=Z(:,noise_start:noise_end)-(meantemp*ones(1,noise_end-noise_start+1));
Znorm=Z;
%Znorm=Znorm/max(max(Znorm));                    % normalize by its maximum value

[max2,dd] = max(Znorm,[],2);  
[fmax,g]=max(max2);
Xmax = dd(g);
XmaxCA = X(Xmax);                               % find the delay corresponding to the maximum scattered power




