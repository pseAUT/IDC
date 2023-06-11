function J=obj(p)
  
  [t,x]=solvemodel(p);

  
  J=x(end,9);

  
end
