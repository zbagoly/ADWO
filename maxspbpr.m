% You are welcome to use these program for your own scientific purposes.  
% If you do so we would appreciate it if you include a reference to the paper
% announcing it: 
% Searching for electromagnetic counterpart of LIGO gravitational waves in the Fermi GBM data with ADWO 
% Z. Bagoly et al., submitted to A&A Letters, http://arxiv.org/abs/1603.06611
%
% maxspbpr returns the Signal Peak to Background Peak Ratio for a given inputvec 
%
% Using ATLAS speeds up the code!

function ret=maxspbpr(inputvec)
        global signal_window; global background_window; global smooth_data;
        inpvec=abs(inputvec); ewm=inpvec(1:6); dw=inpvec(7:20); dw=dw/sum(dw); ewm=ewm/sum(ewm);
        workx=smooth_data*single(reshape(ewm'*dw,84,1));
        ret=-max(workx(signal_window))/max(workx(background_window));
end
