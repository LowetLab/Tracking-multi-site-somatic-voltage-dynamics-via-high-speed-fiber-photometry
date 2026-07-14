% eegfilt() -  (high|low|band)-iass filter data using two-way least-squares 
%              FIR filtering. Multiple data channels and epochs supported.
%              Requires the MATLAB Signal Processing Toolbox.
% Usage:
%  >> [smoothdata] = eegfilt(data,srate,locutoff,hicutoff);
%  >> [smoothdata,filtwts] = eegfilt(data,srate,locutoff,hicutoff, ...
%                                             epochframes,filtorder);
% Inputs:
%   data        = (channels,frames*epochs) data to filter
%   srate       = data sampling rate (Hz)
%   locutoff    = low-edge frequency in pass band (Hz)  {0 -> lowpass}
%   hicutoff    = high-edge frequency in pass band (Hz) {0 -> highpass}
%   epochframes = frames per epoch (filter each epoch separately {def/0: data is 1 epoch}
%   filtorder   = length of the filter in points {default 3*fix(srate/locutoff)}
%   revfilt     = [0|1] reverse filter (i.e. bandpass filter to notch filter). {0}
%
% Outputs:
%    smoothdata = smoothed data
%    filtwts    = filter coefficients [smoothdata <- filtfilt(filtwts,1,data)]
%
% See also: firls(), filtfilt()
%
% Author: Scott Makeig, Arnaud Delorme, SCCN/INC/UCSD (EEGLab Toolbox)

function [smoothdata,filtwts] = eegfilt(data,srate,locutoff,hicutoff,epochframes,filtorder, revfilt)

if nargin<4
    fprintf('');
    help eegfilt
    return
end

if ~exist('firls')
   error('*** eegfilt() requires the signal processing toolbox. ***');
end

[chans frames] = size(data);
if chans > 1 & frames == 1,
    help eegfilt
    error('input data should be a row vector.');
end
nyq            = srate*0.5;  % Nyquist frequency
MINFREQ = 0;

minfac         = 3;    % this many (lo)cutoff-freq cycles in filter 
min_filtorder  = 15;   % minimum filter length
trans          = 0.15; % fractional width of transition zones

if locutoff>0 & hicutoff > 0 & locutoff > hicutoff,
    error('locutoff > hicutoff ???\n');
end
if locutoff < 0 | hicutoff < 0,
   error('locutoff | hicutoff < 0 ???\n');
end

if locutoff>nyq,
    error('Low cutoff frequency cannot be > srate/2');
end

if hicutoff>nyq
    error('High cutoff frequency cannot be > srate/2');
end

if nargin<6
   filtorder = 0;
end
if nargin<7
   revfilt = 0;
end

if isempty(filtorder) | filtorder==0,
   if locutoff>0,
     filtorder = minfac*fix(srate/locutoff);
   elseif hicutoff>0,
     filtorder = minfac*fix(srate/hicutoff);
   end
     
   if filtorder < min_filtorder
        filtorder = min_filtorder;
    end
end

if nargin<5
	epochframes = 0;
end
if epochframes ==0,
    epochframes = frames;    % default
end
epochs = fix(frames/epochframes);
if epochs*epochframes ~= frames,
    error('epochframes does not divide frames.\n');
end

if filtorder*3 > epochframes,   % Matlab filtfilt() restriction
    fprintf('eegfilt(): filter order is %d. ',filtorder);
    error('epochframes must be at least 3 times the filtorder.');
end
if (1+trans)*hicutoff/nyq > 1
	error('high cutoff frequency too close to Nyquist frequency');
end;

if locutoff > 0 & hicutoff > 0,    % bandpass filter
    if revfilt
    else 
    end; 

    f=[MINFREQ (1-trans)*locutoff/nyq locutoff/nyq hicutoff/nyq (1+trans)*hicutoff/nyq 1]; 
    m=[0       0                      1            1            0                      0]; 
elseif locutoff > 0                % highpass filter
 if locutoff/nyq < MINFREQ
    error(sprintf('eegfilt() - highpass cutoff freq must be > %g Hz\n\n',MINFREQ*nyq));
 end
 f=[MINFREQ locutoff*(1-trans)/nyq locutoff/nyq 1]; 
 m=[   0             0                   1      1];
elseif hicutoff > 0                %  lowpass filter
 if hicutoff/nyq < MINFREQ
    error(sprintf('eegfilt() - lowpass cutoff freq must be > %g Hz',MINFREQ*nyq));
 end
 f=[MINFREQ hicutoff/nyq hicutoff*(1+trans)/nyq 1]; 
 m=[     1           1              0                 0];
else
    error('You must provide a non-0 low or high cut-off frequency');
end
if revfilt
    m = ~m;
end;

filtwts = firls(filtorder,f,m); % get FIR filter coefficients

smoothdata = zeros(chans,frames);
for e = 1:epochs                % filter each epoch, channel 
    for c=1:chans
      smoothdata(c,(e-1)*epochframes+1:e*epochframes) ...
        = filtfilt(filtwts,1,data(c,(e-1)*epochframes+1:e*epochframes));
    end
end
