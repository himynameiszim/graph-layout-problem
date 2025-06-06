% pcgWithSSORSolver.m
% Function to solve a linear system Ax=b using the Preconditioned
% Conjugate Gradients (PCG) method with an SSOR preconditioner.

function [x, x_history, res_norm_history, iterations_completed, flag] = ssorSolve(A, b, x0, tol, max_iterations, omega)

    flag = 0; 
    x = x0;
    r = b - A*x;
    
    % solve P_SSOR * z = r
    z = getSSOR(A, r, omega);
    p = z; % first search direction
    rz_old = r'*z; % analogous to residual r'*r of standard CG

    x_history = {x0};
    res_norm_history = [norm(r)];
    iterations_completed = 0;

    for k = 1:max_iterations
        if res_norm_history(end) < tol
            flag = 1;
            return;
        end
        
        % normal CG
        Ap = A*p;
        alpha = rz_old / (p'*Ap);        
        x = x + alpha*p;
        r_new = r - alpha*Ap; % Store new residual temporarily
        
        x_history{end+1} = x;
        res_norm_history(end+1) = norm(r_new);
        iterations_completed = k;

        % apply preconditioner
        z_new = getSSOR(A, r_new, omega);
        rz_new = r_new'*z_new;
        beta = rz_new / rz_old;
        p = z_new + beta*p; % new search direction
        
        r = r_new;
        rz_old = rz_new;
    end
end
