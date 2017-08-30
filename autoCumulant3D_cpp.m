function [ cum_img ] = autoCumulant3D_cpp( imgstack, order, time_lags  )
% USAGE: [ cum_img ] = autoCumulant3D_cpp( imgstack, order, time_lags  )
% Computes the cumulant 'cum' of order 'order' (max 4) for the 4D (x,y,z,
% time) image stack 'imgstack'  with time lags 'time_lags'. The required 
% number of specified time lags is #(order-1).
%
% Output
% cum_img(xdim, ydim, zdim) - the cumulant for every pixel
%
%
% Author: Simon Christoph Stein
% Date: 2017

end

