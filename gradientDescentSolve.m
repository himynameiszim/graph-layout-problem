% gradientDescentSolve.m
% Function to solve a linear system Ax=b using the Gradient Descent method.

function [x, x_old, residual_norm_old, total_iterations, flag] = gradientDescentSolve(A, b, x0, eta, tol, max_iterations)
    
    flag = 0;
    x = x0; % intial guess
    x_old = {x0}; % keep track of x each iteration
    residual_norm_old = [norm(b - A*x)]; % keep track of residual(s) each iteration
    total_iterations = 0; % keep track of number of iterations until convergence (or halt if reach max_iterations)

    for k=1:max_iterations
        if residual_norm_old(end) < tol
            flag = 1;
            break;
        end
        
        gradient = 2*(A*x - b);
        x = x - eta*gradient;
        
        x_old{end+1} = x;
        residual_norm_old(end+1) = norm(b - A*x);
        total_iterations = k;
    end
end