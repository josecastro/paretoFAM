%%% Implementation of pareto FSFAM ARTMAP in Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% The input and output

% TODO
%   1. Comment Code
%   2. Create Test Cases for quickstart, quicksortRest
%   3. Create Test Cases for pareto front generation
%


load -binary 'datos.bin'
% loads num_pats, a, bmat

% NEXT LINE COMMENTED
% bmat is only one value of categories
% bmat = [ bmat; 1-bmat ];        % Change to [1 0], [0 1] form

ac = [a; ones(size(a))-a];   % The complement-coded form of input a

start = time();

%%%%%%%%%%%%%%%%% If we are testing, make a grid of test inputs
test_mode = 0;    % Zero means we're training, 1 means testing

if test_mode==1,
        grain = 100;
        grx = linspace(xcent-sq/2,xcent+sq/2,grain)';
        gry = linspace(ycent-sq/2,ycent+sq/2,grain)';
        testpats = zeros(grain^2,2);    % Initialise
        for i=1:grain,
                testpats(1+grain*(i-1):i*grain,1) = grx;
                testpats(1+grain*(i-1):i*grain,2) = ...
                                    ones(grain,1)*gry(i);
        end;
        testpats = [ testpats' ; 1-testpats' ];
        test_out = zeros(size(testpats,2),1);   
                             % Initialise output vector
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters

alpha = 0.001;  % "Choice" parameter > 0. Set small for the
                % conservative limit (Fuzzy AM paper, Sect.3)
beta = 1;       % Learning rate. Set to 1 for fast learning
rho_bar = 0;    % Baseline vigilance for ARTa, in range [0,1]
M = size(a,1);  % Number of input components. Derived from data
                % NB: Total input size = 2M (due to complement)
N = 1;         % Number of available coding nodes
                % We start with some resonably large number
                % then, if we need to, can add more uncommitted
L = size(bmat,1);       % Number of output nodes. Derived from data
epsilon = 0.001;        % Fab mismatch raises ARTa vigilance to this 
                        % much above what is needed to reset ARTa
num_patterns = size(a,2);   % Number of input patterns

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Set up weights

w = [ones(2*M,1); -1];  % Initial weights in ARTa. All set to 1,
                        % classification is none = -1 uncommitted node
                        % Row-i, col-j entry = weight from input node i
                        % to F2 coding node j
class = 2*M+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Main loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

for i=1:num_patterns    % Go through patterns one by one 
                        % Note: could shuffle order using randperm 

        A = ac(:,i);    % Present input is i-th column of ac

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
        %%%%%%%%%% Find the winning, matching ARTa node
                
        N = size(w,2);   
        % Count how many F2a nodes we have

        A_for_each_F2_node = A * ones(1,N);
        % Matrix containing a copy of A for 
        % each F2 node. Useful for Matlab

        A_AND_w = min(A_for_each_F2_node,w(1:2*M,:));  
        % Fuzzy AND = min

        S = sum(A_AND_w);       
        % Row vector of signals to F2 nodes

        rho_vect = S / sum(A);
        candidate_indexes = find(rho_vect >= rho_bar);
        rho_vect = rho_vect(:,candidate_indexes);
        T_vect = S(:,candidate_indexes) ./ (alpha+sum(w(1:2*M,candidate_indexes)));
        
        TSort = [T_vect' (1:length(T_vect))'];
        [TSort stack] = quickstart(TSort);
        index = candidate_indexes(TSort(1,2));
        
        if w(class,index) ~= -1 && w(class,index) ~= bmat(i)
          TSort = quicksortRest(TSort, stack);
          T_ind = TSort(:,2)';
          
          [_,rho_ind] = sort(rho_vect,'descend');

          rho_vect = candidate_indexes(rho_ind)';
          T_vect   = candidate_indexes(T_ind)';
          rhoFirst = rho_ind(1);
          TFirst   = T_ind(1);
          rhoMat = [candidate_indexes(rho_ind)' [rho_ind(2:length(rho_ind))' ; 0] [0; rho_ind(1:length(rho_ind)-1)']];
          TMat   = [candidate_indexes(T_ind)'   [T_ind(2:length(T_ind))' ; 0]     [0; T_ind(1:length(T_ind)-1)']];
          rhoMat(rho_ind,:) = rhoMat;
          TMat(T_ind,:) = TMat;
          
          prox  = 2;

          index = TMat(TFirst,1);
          while index ~= 0 && w(class,index) ~= -1 && w(class,index) ~= bmat(i)
            [N TMat rhoMat] = elimineResto(rhoMat(TFirst,prox), TMat, rhoMat);
            [_ rhoMat] = elimine(TFirst, rhoMat);
            [TFirst, TMat] = elimine(TFirst, TMat);
            if TFirst ~= 0
              index = TMat(TFirst,1);
            else
              index = 0;
            end
          end
        end
        
        if (w(class,index) == -1) % uncommited node
          w = [w [A; bmat(i)]];
        elseif (w(class,index) == bmat(i))
          w(1:2*M,index) = (1-beta)*w(1:2*M,index) + beta*min(w(1:2*M,index), A);
        else
          'Error'
        end
        
end

[_ sort_order] = sort(w(1,:));
w = w(:,sort_order);
disp('iteration: '), disp(i);
%disp(w);
time() - start
disp('finished');


