function yi = fastInterpTable2D(X1, X2, Y, x1i, x2i)
   if isempty(x1i) || isempty(x2i)
       yi = [];
       return
   end
   size_x = size(x1i);
   if size_x(1) == 1 || size_x(2) == 1
       [x1i, x2i] = ndgrid(x1i, x2i); 
   end

   F = griddedInterpolant(X1, X2, Y, 'linear', 'linear');
   yi = F(x1i, x2i);
   disp('fastInterp yi:')
   disp(size(yi))
end