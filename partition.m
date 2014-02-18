%
% Matlab quicksort utility function
%   Each data is a row in a matrix
%   asuming key is in first column of matrix
%   the rest of the row corresponds to data
%

function [storeIndex array] = partition(array, left, right)
   pivotIndex = idivide(left+right,2,'fix');
   pivotValue = array(pivotIndex,1);

   temp = array(pivotIndex,:);
   array(pivotIndex,:) = array(right,:);
   array(right,:) = temp;
   
   storeIndex = left;
   for i = left:(right-1)
     if array(i,1) >= pivotValue
       temp = array(i,:);
       array(i,:) = array(storeIndex,:);
       array(storeIndex,:) = temp;
       
       storeIndex = storeIndex+1;
     end
   end
   
   temp = array(storeIndex,:);
   array(storeIndex,:) = array(right,:);
   array(right,:) = temp;
endfunction
   
   
