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



graph

;%% Min z = 3x1 - 5x2 
% where s.t. are 
% x1 + x2 <= 6 
% 2x1 - x2 >= 9 
% x1 + 2x2 <= 6
% x1 >= 0
% x2 >= 0

clc ;
clear;

%% Phase 1, input parameter
C = [3 -5];
Z = @(x1, x2)3*x1 - 5*x2;

% Define the constraints
A = [1 1; 2 -1; 1 2];
b = [6; 9; 6];

% Defining st's as functions as well 
C1 = @(x1, x2) x1 + x2 - 6;
C2 = @(x1, x2) 2*x1 - x2 - 9;
C3 = @(x1, x2) x1 + 2*x2 - 6;

% finding number of constraints 
[m, n] = size(A);

%% Phase 2 Plotting
x1 = 0:max(b./A(:, 1)); % for x axis so that we dont have to plot a huge graph, only what is required
% Known as feasable region
for i = 1:m
    x2 = (b(i) - A(i, 1) * x1) / A(i, 2);

    plot(x1, x2)
    hold on
end 

%% Phase 3 Find intersection & corner points
% adding st's for x1 and x2 to be positive 
A = [A; eye(n)];
b = [b; zeros(2, 1)]; % adding constraint for first quadrant

m = size(A, 1); % updating
pt = [];
for i = 1:m
    for j = i + 1: m 
        aa = [A(i, :); A(j, :)];
        bb = [b(i); b(j)];
        if (det(aa))
            X = aa \ bb;
            if(X >= 0)
                pt = [pt X];
            end
        end
    end
end

pt = unique(pt', "rows")';
disp(pt)

%% Phase 4 Find Feasable points 
FP = [];
z = [];
for i = 1:size(pt, 2) % 2 gives columns, i.e. loop in number of points
    pt1 = pt(1, i); % x1 of ith point
    pt2 = pt(2, i); % x2 of ith point
    if(C1(pt1, pt2) <= 0 && C2(pt1, pt2) >= 0 && C3(pt1, pt2) <= 0)
        FP = [FP pt(:, i)]; % Store feasible points
        plot(pt1, pt2, '*r', 'MarkerSize', 30);
        cost = Z(pt1, pt2);
        z = [z cost];
    end
end

disp(FP)
disp(z)

%% Phase 5 finding optimal solution and value
[optimal_val index] = min(z)
optimal_sol = FP(:, index)
% Display the optimal solution and its corresponding value
fprintf('Optimal Solution: x1 = %.2f, x2 = %.2f\n', optimal_sol(1), optimal_sol(2));
fprintf('Optimal Value: Z = %.2f\n', optimal_val);










