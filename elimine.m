%
% Matlab/Octave function to eliminate a value from a
% cursor list, implemented as an array
% size(vect) = (N,3)
%    vect(K,1) = value
%    vect(K,2) = next entry
%    vect(K,3) = prev entry
%    0 means end of list
%

function [Val vect] = elimine(pos, vect)
   id   = 1;
   prox = 2;
   ant  = 3;
   
   if (vect(pos,ant) ~= 0)
     vect(vect(pos,ant),prox) = vect(pos,prox);
   end
   if (vect(pos,prox) ~= 0)
     vect(vect(pos,prox),ant) = vect(pos,ant);
   end
   Val = vect(pos,prox);
end
