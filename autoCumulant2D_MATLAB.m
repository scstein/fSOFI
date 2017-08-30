function cum = autoCumulant2D_MATLAB(imgstack,order,time_lags)
% USAGE: [ cum_img ] = autoCumulant2D_cpp( imgstack, order, time_lags  )
% Computes the cumulant 'cum' of order 'order' (max 4) for the 3D (x,y,time)
% image stack 'imgstack'  with time lags 'time_lags'. The required 
% number of specified time lags is #(order-1).
%
% Output
% cum(xdim, ydim) - the cumulant for every pixel
%
%
% Author: Simon Christoph Stein
% Date: 2017

  if numel(time_lags) ~= order-1
    error('Need #(order-1) time lags.');
  end
  
  max_tlag = max(time_lags);
  nFrames = size(imgstack,3)-max_tlag;
  
  if max_tlag >= size(imgstack,3)
      error('Max. time lag needs to be smaller than the number of frames in the movie!');
  end
  
  % subtract avg;
  avg = mean(imgstack,3);
  imgstack = imgstack - repmat(avg,[1,1,size(imgstack,3)]);
  

switch order
    case 1
        cum = avg;
    case 2 
        cum = mean(imgstack(:,:,1:(1+nFrames-1)).*imgstack(:,:,1+time_lags(1):(1+time_lags(1)+nFrames-1)),3);
    case 3        
        cum = mean(imgstack(:,:,1:(1+nFrames-1)).*imgstack(:,:,1+time_lags(1):(1+time_lags(1)+nFrames-1)).*imgstack(:,:,1+time_lags(2):(1+time_lags(2)+nFrames-1)),3);
    case 4
        cum = -mean(imgstack(:,:,1:(1+nFrames-1)).*imgstack(:,:,1+time_lags(3):(1+time_lags(3)+nFrames-1)),3).*mean(imgstack(:,:,1+time_lags(1):(1+time_lags(1)+nFrames-1)).*imgstack(:,:,1+time_lags(2):(1+time_lags(2)+nFrames-1)),3) ...
              -mean(imgstack(:,:,1:(1+nFrames-1)).*imgstack(:,:,1+time_lags(2):(1+time_lags(2)+nFrames-1)),3).*mean(imgstack(:,:,1+time_lags(1):(1+time_lags(1)+nFrames-1)).*imgstack(:,:,1+time_lags(3):(1+time_lags(3)+nFrames-1)),3) ...
              -mean(imgstack(:,:,1:(1+nFrames-1)).*imgstack(:,:,1+time_lags(1):(1+time_lags(1)+nFrames-1)),3).*mean(imgstack(:,:,1+time_lags(2):(1+time_lags(2)+nFrames-1)).*imgstack(:,:,1+time_lags(3):(1+time_lags(3)+nFrames-1)),3) ...
              +mean(imgstack(:,:,1:(1+nFrames-1)).*imgstack(:,:,1+time_lags(1):(1+time_lags(1)+nFrames-1)).*imgstack(:,:,1+time_lags(2):(1+time_lags(2)+nFrames-1)).*imgstack(:,:,1+time_lags(3):(1+time_lags(3)+nFrames-1)),3);            
end


return


