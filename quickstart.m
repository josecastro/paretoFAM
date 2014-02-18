%
% Matlab/Octave implementation of lazy quicksort
%   used for optimistic ordering when you search
%   through an ordered list but might need only
%   the first element
%

function [array stack] = quickstart(array)
   low = 1;
   high = size(array,1);
   stack = [[0 -1]];
   while low <= high
     [pivotLocation array] = partition(array, low, high);
     stack = [[pivotLocation+1 high] ; stack];
     high = pivotLocation - 1;
   end
endfunction
