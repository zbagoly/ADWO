% You are welcome to use these program for your own scientific purposes.  
% If you do so we would appreciate it if you include a reference to the paper
% announcing it: 
% Searching for electromagnetic counterpart of LIGO gravitational waves in the Fermi GBM data with ADWO 
% Z. Bagoly et al., submitted to A&A Letters, http://arxiv.org/abs/1603.06611
%
% Using ATLAS speeds up the code!

%INPUT :  smooth_input_data tensor : 	should be provided, it is 1ms resolution smoothed data, 
% 					every time slice is a (channels x detectors) matrix
disp(size(smooth_input_data)) 		% should be (700001, 6, 14);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
polinome_order=6; %%% background fit polynome order 

adwotime=linspace(-200,500,700*1000+1); % the total time of the smoothed data, 1ms step 
NT=size(adwotime,2); 
smooth_data=single(zeros(NT,6*14));  	% single precision is for a speedup in octave

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data reshape start
for detector=1:14 %% all detectors (14)
  for channel=1:6 %% all energy channels (6)
	work_lc=reshape(smooth_input_data(:,channel,detector),NT,1);  
	fit_parameters=polyfit(adwotime,work_lc',polinome_order);   		% polynome fit 
	smooth_data(:,indx)=single(wx'-polyval(fit_parameters,adwotime)); 	% background substracted smooth data
	indx=indx+1;
  end
end
global smooth_data; % should be checked in octave - sometimes the global flag is not set from the command line
% data reshape end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADWO setup start
run_positions=[-192:6:492]; 	% whole 6s windows, covering the interval, centered on 0
window_center=0; signal_window=(find((adwotime<=window_center+3) .* (adwotime>=window_center-3))); nt3=adwotime*0;nt3(signal_window)=1; background_window=find(nt3==0); 
global signal_window; global background_window; 
indx=1; maxspbprout=zeros(size(run_positions,1),1+1+6+14);  % output matrix
% ADWO setup end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADWO loop start
% We invesigate every 6s window in the full adwotime region and optimize the S/B Peak ratio for each window
% The peak postion within the signal_window is a free parameter!
%
for window_center=run_positions
	signal_window=(find((adwotime<=window_center+3) .* (adwotime>=window_center-3))); nt3=adwotime*0;nt3(signal_window)=1; background_window=find(nt3==0); 
	global signal_window; global background_window; global smooth_data; % should be checked in octave - sometimes the global flag is not set from the command line
	edweights=[ones(1,6)/6 ones(1,14)/14]; % optimalization init
	% ADWO : altough fminsearch reports the maxima, restarting iproves it slightly. 
	% Should be optimized by comparing the results in each step (in the next versions)
	% Now we restart around the optimum several times: after it converged, it is fast and cheap to restart again ...
	optweights=abs(edweights); [edweights maxspbpr]=fminsearch(@maxspbpr,optweights,optimset('MaxIter',82500,'TolX', 1e-5)); 
	optweights=abs(edweights); [edweights maxspbpr]=fminsearch(@maxspbpr,optweights,optimset('MaxIter',82500,'TolX',1e-10));
	optweights=abs(edweights); [edweights maxspbpr]=fminsearch(@maxspbpr,optweights,optimset('MaxIter',82500,'TolX',1e-10));
	optweights=abs(edweights); [edweights maxspbpr]=fminsearch(@maxspbpr,optweights,optimset('MaxIter',82500,'TolX',1e-12));
	% ADWO converged, now do the weight normalization
	optweights=abs(edweights); cw=optweights(1:6); dw=optweights(7:20); dw=dw/sum(dw); cw=cw/sum(cw); 
	% fill up the output matrix: the peak is somewhere in the signal_window! 
	maxspbprout(indx,:)=[window_center abs(maxspbpr) cw dw]; 
	indx=indx+1;
end
% ADWO loop end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% detector weights and ADWO plot with a window around T=0 
indx=find(maxspbprout(:,1)==0);
ew=maxspbprout(indx,3:8); dw=maxspbprout(indx,9:22));
workx=smooth_data*single(reshape(ew'*dw,84,1));
disp(ew,dw); 
plot(adwotime,workx,'-');

