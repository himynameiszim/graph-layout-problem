% conjugateGradientSolve.m
% Function to solve a linear system Ax=b using the Conjugate Gradients method.
% This is the standard, computationally efficient implementation.
% A must be symmetric and positive-definite.

function [x, x_history, res_norm_history, total_iterations, flag] = conjugateGradientSolve(A, b, x0, tol, max_iterations)
    
    flag = 0;
    x = x0;
    r = b - A*x;
    p = r; % search direction
    
    rs_old = r'*r; 
    x_history = {x0};
    res_norm_history = [norm(r)];
    total_iterations = 0;
    
    for k = 1:max_iterations
        if res_norm_history(end) < tol
            flag = 1;
            return;
        end        
        Ap = A*p; 
        
        alpha = rs_old / (p'*Ap); % step size
        x = x + alpha*p; % update estimate
        r = r - alpha*Ap; % update residual
        rs_new = r'*r; % norm of residual
        
        % store history
        x_history{end+1} = x;
        res_norm_history(end+1) = sqrt(rs_new);
        total_iterations = k;
        
        beta = rs_new / rs_old;
        p = r + beta*p;
        rs_old = rs_new;
    end
end
