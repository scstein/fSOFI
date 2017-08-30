function [sof, win_sof] = FSOFI_Analysis3D(im, ncum, win_size, time_lags, itp_factor, mirrorMode)
% Standard usage: sof = SOFIAnalysis(im, order, win_size, max_itp_factor)
%
% Computes SOFI images up to the specified order from a movie 'im' with fluctuating emitters.
%
% Input
% im - 4D (x,y,z, time) image stack
% ncum - SOFI order to compute
% win_size - Size of data window used to compute the output | default: 100
% time_lags - vector with #(order-1) time lags OR 0, which sets all lags to
%             zero. | default: [1,2,3,...]
% itp_factor - Magnification factor for subpixel generation for each
%              dimension [ipX,ipY,ipZ]. If a single number is given, the same
%              factor is used for all dimensions. The output is of size
%              itp_fac.*size(img)| default: 1
% mirrorMode - Specifies wether to use periodic padding of the input data before
%           fourier interpolation (example see fSOFI publication).
%           The padding prevents artifacts from non-periodic borders and is
%           essential if only a low number of pixels is available along a
%           specific dimension. Possible values: 'none','lateral','axial','both'
%           | default: 'none'
%
% Output
% sof(xdim, ydim,zdim) - SOFI result image
% win_sof(xdim,ydim,zdim,nr_windows) - SOFI result for every window
%
% Author: Simon Christoph Stein
% E-Mail: scstein@phys.uni-goettingen.de
% Date: 2017

if nargin <3 || isempty(win_size)
    win_size = min( size(im,4), 100);
else
    win_size = min( size(im,4), win_size);
end

if nargin <4 || isempty(time_lags)
    time_lags = 1:(ncum-1); % default time lags 1,2,3,...
else
    if  time_lags(1) == 0 && numel(time_lags) == 1 % Special case for 0 time lag
        time_lags = repmat(time_lags, 1, ncum-1);
    end
end

NrDIMS = 3;
if nargin <5 || isempty(itp_factor)
    itp_factor = 1;
end
if ~( numel(itp_factor) == 1 || numel(itp_factor) == NrDIMS)
    error('%i interpolation factors specified. Give either one for all dimension or one per dimension!', numel(itp_factor))
end
if numel(itp_factor) == 1
    itp_factor = repmat(itp_factor,[1,NrDIMS]);
end

if nargin <6 || isempty(mirrorMode)
    mirrorMode = 'none';
end

if ncum>4
    error('SOFIAnalysis:argChk', 'Higher cumulants than 4th order are not supported.')
end

% im = double(im); % Image must be double or single

full_time = tic;

% Memory allocation
win_nr = floor(size(im,4)/win_size); % Number of windows
ip_im = zeros( itp_factor(1)*size(im,1), itp_factor(2)*size(im,2), itp_factor(3)*size(im,3), win_size); % Storage for interpolated images
win_sof = zeros(size(ip_im,1),size(ip_im,2), size(ip_im,3), win_nr); % Result for each individual window
sof = zeros(size(ip_im,1),size(ip_im,2), size(ip_im,3));             % Final SOFI result


msgAccumulator = ''; % Needed for rewindable command line printing (rewPrintf subfunction)
try
    for iWin = 1:win_nr
        rewindMessages(); % Remove cached messages from command line
        rewPrintf('Window %i/%i \n', iWin, win_nr);
        
        min_frame = (iWin-1)*win_size + 1;
        max_frame = iWin*win_size;
        
        % Fourier interpolation if requested
        if(max(itp_factor) > 1)
            rewPrintf('  Fourier interpolation.. \n');
            
            cnt = 1;
            % Here we fill the container ip_im with interpolated images from the current window.
            for iFrame = min_frame:max_frame
                ip_im(:,:,:,cnt) = fourierInterpolation( double(im(:,:,:,iFrame)), itp_factor, mirrorMode);
                cnt = cnt+1;
            end
            rewPrintf('done\n');
            
            %--- Core SOFI Computation ---
            rewPrintf('  Computing the cumulant.. ');
            win_sof(:,:,:,iWin) = autoCumulant3D_cpp(ip_im,ncum, time_lags); % C++ (recommended)
            %             win_sof(:,:,:,iWin) = autoCumulant3D_cpp(ip_im,ncum, time_lags); % MATLAB
            sof = sof + win_sof(:,:,:,iWin);
            rewPrintf('done\n');
        else
            %--- Core SOFI Computation ---
            rewPrintf('  Computing the cumulant.. ');
            win_sof(:,:,:,iWin) = autoCumulant3D_cpp(double(im(:,:,:,min_frame:max_frame)),ncum, time_lags); % C++ (recommended)
            %             win_sof(:,:,:,iWin) = autoCumulant3D_cpp(double(im(:,:,:,min_frame:max_frame)),ncum, time_lags); % MATLAB
            sof = sof + win_sof(:,:,:,iWin);
            rewPrintf('done\n');
        end
    end
    sof = sof/win_nr;
catch Merr
    warning('Error during FSOFI_Analysis3D. Did you install the ''Microsoft Visual Studio 2012 Redistributable Version 4'' for the MEX functions to work? If not, please install (recommended) or uncomment the MATLAB version of the cumulant calculation inside the code (not recommended).');
    throw(Merr);
end

toc(full_time);

return


    function rewPrintf(msg, varargin)
        % Rewindable message printing: Print msg and cache it.
        % Usage is analogous to sprintf.
        msg = sprintf(msg, varargin{:});
        msgAccumulator = [msgAccumulator, msg];
        fprintf(msg);
    end

    function rewindMessages()
        % Remove cached messages from command line, reset cache
        reverseStr = repmat(sprintf('\b'), 1, length(msgAccumulator));
        fprintf([reverseStr]);
        
        msgAccumulator = '';
    end
end

