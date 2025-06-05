clear; clc; close all;

% 1. graph setup
total_vertices = 5; 
edges = [1, 2;
         2, 3;
         3, 1;
         1, 4;
         2, 5;
         4, 5];
V_idx = (1:total_vertices);
V0_idx = [3; 4]; % specified vertices
U_idx = setdiff(V_idx, V0_idx); % unknown vertices
m = length(U_idx);
V_position = [
    2.0, 3.0;  
    2.5, 0.5;  
    0.0, 0.0; 
    4.0, 0.0; 
    1.5, 3.5  
]; % initial position of all vertices

% 2. setting up linear systems
[L_mat, degrees_mat, adjacency_mat] = buildGraphLaplacian(total_vertices, edges, U_idx);

% init x_0, y_0
x_init_guess = V_position(U_idx, 1);
y_init_guess = V_position(U_idx, 2);

% construct b_x and b_y
bx = zeros(m, 1);
by = zeros(m, 1);
for k=1:m
    neighbours_of_k = adjacency_mat{U_idx(k)};
    fixed_neighbour_x = 0;
    fixed_neighbour_y = 0;
    for neighbour_idx=1:length(neighbours_of_k)
        neighbour = neighbours_of_k(neighbour_idx);
        if ismember(neighbour, V0_idx) %if neighbour of unknown vertex is fixed
            fixed_pos_neighbour = V_position(neighbour, :);
            fixed_neighbour_x = fixed_neighbour_x + fixed_pos_neighbour(1);
            fixed_neighbour_y = fixed_neighbour_y + fixed_pos_neighbour(2);
        end
    end
    bx(k) = fixed_neighbour_x;
    by(k) = fixed_neighbour_y;
end

% 3. set up gradient descent
tol = 1e-6;
max_iteration = 500;
eta = 0.05;

% 4. gradient descent solve
[x_gd_final, x_gd_old, res_norm_x_old, x_num_iteration, x_flag] = gradientDescentSolve(L_mat, bx, x_init_guess, eta, tol, max_iteration);
[y_gd_final, y_gd_old, res_norm_y_old, y_num_iteration, y_flag] = gradientDescentSolve(L_mat, by, y_init_guess, eta, tol, max_iteration);
if x_flag * y_flag == 1
    fprintf('Gradient descent converged in %d iteration(s)', max(x_num_iteration, y_num_iteration));
end
if x_flag * y_flag == 0
    fprintf('Gradient descent diverged after %d iteration(s)', max_iteration);
end









