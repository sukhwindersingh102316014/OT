bfs 

clc; 
clear all;
clear figure;

C = [1 2 0 0];

A = [-1 1 1 0;
      1 1 0 1];   % Constraint coefficients

b = [1; 2];       % RHS

z = @(X) C*X;

m = size(A,1);    % number of constraints
n = size(A,2);    % number of variables

%% Phase 2: Basic Solutions and BFS
basicsol = [];
bfsol = [];

% >>> ADDED MATRICES <<<
degenerateBFS = [];
generateBFS   = [];

pair = nchoosek(1:n, m);

for i = 1:size(pair,1)

    basicvar_index = pair(i,:);
    y = zeros(n,1);

    % Solve for basic variables
    X = A(:, basicvar_index) \ b;
    y(basicvar_index) = X;

    basicsol = [basicsol y];

    % -------- Feasibility Check --------
    if all(X >= 0)

        bfsol = [bfsol y];

        % -------- Degenerate / Non-degenerate Check --------
        if sum(y > 0) < m
            degenerateBFS = [degenerateBFS y];
        else
            generateBFS = [generateBFS y];
        end
    end
end

disp('All Basic Solutions:');
disp(basicsol)

disp('Basic Feasible Solutions:');
disp(bfsol)

disp('Degenerate BFS:');
disp(degenerateBFS)

disp('Non-Degenerate BFS:');
disp(generateBFS)

%% Phase 3: Optimal Solution
cost = z(bfsol);
[opt_val, index] = max(cost);
optsol = bfsol(:, index);

disp('Optimal Solution:');
disp(optsol)

disp('Optimal Value:');
disp(opt_val)












