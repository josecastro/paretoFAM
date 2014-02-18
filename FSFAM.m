%%% Script Implementation of FSFAM ARTMAP in Matlab/Octave

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% The input and output

% TODO
%   1. Comment Code
%   2. Create Test Cases for FSFAM
%

load -binary 'datos.bin'

start = time();

ac = [a; ones(size(a))-a];   % The complement-coded form of input a

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

class_idx = 2*M+1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Set up weights

w = [ones(2*M,1); -1];  % Initial weights in ARTa. All set to 1,
                        % classification is none = -1 uncommitted node
                        % Row-i, col-j entry = weight from input node i
                        % to F2 coding node j

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

        %%%%%%%%%%%%%%%%% Initialise
        rho = rho_bar;      
        % We start off with ARTa vigilance at baseline
        resonant = 0;         
        % We're not resonating in the ARTa module yet

        rho_vect = S / sum(A);
        ind_pass = find(rho_vect >= rho);
        rho_pass = rho_vect(ind_pass);
        
        NMatch = 0;
        while (~resonant)
          T = rho_pass * sum(A) ./ (alpha+sum(w(1:2*M,ind_pass)));
          [ Tmax, J ] = max(T);
          rho_now = rho_pass(J);
          J = ind_pass(J);
          if (w(class_idx,J) == -1) % uncommited node
            w = [w, [A; bmat(i)]];
            resonant = 1;
          elseif (w(class_idx,J) == bmat(i))
            w(1:2*M,J) = (1-beta)*w(1:2*M,J) + beta*min(w(1:2*M,J), A);
            resonant = 1;
          else
            % matchtracking
            NMatch   = NMatch+1;
            rho      = rho_now + epsilon;
            rho_vect = S / sum(A);
            ind_pass = find(rho_vect >= rho);
            rho_pass = rho_vect(ind_pass);
          end
        end
end

[_ sort_order] = sort(w(1,:));
w = w(:,sort_order);
disp('iteration: '), disp(i);
%disp(w);
time() - start
disp('finished.');
