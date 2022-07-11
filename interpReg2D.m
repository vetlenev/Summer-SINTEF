function yi = interpReg2D(T, x1i, x2i, reginx)
   if isempty(x1i) || isempty(x2i)
      yi     = [];     
      return;
   end
   nreg = numel(reginx);

   if nreg > 0 
      if ischar(reginx{1}) && strcmp(reginx{1}, ':')
          % Special case denoting entire domain in single region.
          yi = fastInterpTable2D(T{1}{1}, T{1}{2}, T{1}{3}, x1i, x2i);
      else
          % Iterate over region indices and assign values
          yi = zeros(numel(x1i), numel(x2i));
          for k = 1:nreg
             subs = reginx{k};
             if ~isempty(subs)                 
                disp('interpReg yi:')
                disp(size(yi))
                disp(size(yi(subs, subs)))
                yi(subs, subs) = fastInterpTable2D(T{k}{1}, T{k}{2}, T{k}{3}, x1i(subs), x2i(subs));
             end
          end
      end
   end           
end