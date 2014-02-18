%
% Matlab/Octave lazy quicksort implementation
%  uses a stack (matrix) with 
%


function vect = quicksortRest(vect, stck)
   while size(stck,1) > 0
   
     % pop values from the stack
     low = stck(1,1);
     high = stck(1,2);
     stck = stck(2:size(stck,1),:);
     
     while low <= high
       [pivotLocation vect] = partition(vect, low, high);
       
       % push values into the stack
       stck = [[pivotLocation+1 high] ; stck];
       
       high = pivotLocation -1;
     end
   end
endfunction
