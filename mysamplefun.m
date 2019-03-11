function  vec=mysamplefun(m, n)
  % sample n numbers from 1:m
  vec = ceil(rand(1, n) * m);
  return 
  
  while (1)
    x = ceil(rand(1, n*2) * m);
    y = unique(x,'stable'); 
    %y = unique(x);
    if (length(y) >= n)
       break
    end
  end
  vec = y(1:n);
