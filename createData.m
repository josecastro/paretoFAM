%
% Circle in the square data generator
% reference pending ...
%

function dummy = createData(num_pats)
   sq = 1;                         % Size of square
   r = sq/sqrt(2*pi);              % Radius of circle so it's half area of square
   xcent = 0.5;
   ycent = 0.5;                    % Centre of circle
   a = [ xcent*ones(1,num_pats); ycent*ones(1,num_pats)] + sq*(0.5-rand(2,num_pats));              % The x,y coords
   bmat = ((a(1,:)-xcent).^2 + (a(2,:)-ycent).^2) > r^2;
   save -binary 'datos.bin' num_pats a bmat;
endfunction
