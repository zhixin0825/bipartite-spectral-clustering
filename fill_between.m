function   fill_between( x, y1, y2, varargin )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
   x = x(:);
   y1 = y1(:);
   y2 = y2(:);
   
   fill([x; flipud(x)], [ y1; flipud(y2) ], varargin{:}) 
end

