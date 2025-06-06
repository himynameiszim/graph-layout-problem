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

% 3. set up for solve
% init x_0, y_0
x_init_guess = V_position(U_idx, 1);
y_init_guess = V_position(U_idx, 2);
tol = 1e-8;
max_iteration = 500;
eta = 0.05;

% 4. gradient descent solve
[x_gd_final, x_gd_old, resNorm_x_gd_old, x_gd_numIter, x_gd_flag] = gradientDescentSolve(L_mat, bx, x_init_guess, eta, tol, max_iteration);
[y_gd_final, y_gd_old, resNorm_y_gd_old, y_gd_numIter, y_gd_flag] = gradientDescentSolve(L_mat, by, y_init_guess, eta, tol, max_iteration);
if x_gd_flag && y_gd_flag
    fprintf('Gradient descent converged to \n x = [ ');
    fprintf('%.8f ', x_gd_final);
    fprintf(']\n and \n y =[ ');
    fprintf('%.8f ', y_gd_final);
    fprintf(']\n in %d iteration(s) \n\n', max(x_gd_numIter, y_gd_numIter));
    else 
    fprintf('Gradient descent diverged after %d iteration(s) \n\n', max_iteration);
end

% 5. conjugate gradients solve
[x_cg_final, x_cg_old, resNorm_x_cg_old, x_cg_numIter, x_cg_flag] = conjugateGradientSolve(L_mat, bx, x_init_guess, tol, max_iteration);
[y_cg_final, y_cg_old, resNorm_y_cg_old, y_cg_numIter, y_cg_flag] = conjugateGradientSolve(L_mat, by, y_init_guess, tol, max_iteration);
if x_cg_flag && y_cg_flag
    fprintf('Conjugate gradients converged to \n x = [ ');
    fprintf('%.8f ', x_cg_final);
    fprintf(']\n and \n y =[ ');
    fprintf('%.8f ', y_cg_final);
    fprintf(']\n in %d iteration(s) \n\n', max(x_cg_numIter, y_cg_numIter));
    else 
    fprintf('Conjugate gradients diverged after %d iteration(s) \n\n', max_iteration);
end

% 6. SSOR-preconditioned conjugate gradients solve 
omega = 0.5;
[x_pcg_final, x_pcg_old, resNorm_x_pcg_old, x_pcg_numIter, x_pcg_flag] = ssorSolve(L_mat, bx, x_init_guess, tol, max_iteration, omega);
[y_pcg_final, y_pcg_old, resNorm_y_pcg_old, y_pcg_numIter, y_pcg_flag] = ssorSolve(L_mat, by, y_init_guess, tol, max_iteration, omega);
if x_pcg_flag && y_pcg_flag
    fprintf('SSOR preconditioned conjugate gradients converged to \n x = [ ');
    fprintf('%.8f ', x_pcg_final);
    fprintf(']\n and \n y =[ ');
    fprintf('%.8f ', y_pcg_final);
    fprintf(']\n in %d iteration(s) \n\n', max(x_pcg_numIter, y_pcg_numIter));
    else 
    fprintf('Conjugate gradients diverged after %d iteration(s) \n\n', max_iteration);
end




