% getSSORres.m
% function to solve P_SSOR * z = r where P_SSOR is the preconditioner
% matrix
% P_SSOR = (1/(omega*(2-omega))) * (D + omega*L) * inv(D) * (D + omega*L')

function z = getSSOR(A, r, omega)
    % diagonal (D) and strictly lower triangular (L)
    D = diag(diag(A));
    L_slt = tril(A, -1); % Strictly lower triangular part

    % solve P_SSOR * z = r.   
    r_prime = omega * (2 - omega) * r;
    w = (D + omega * L_slt) \ r_prime;
    v = D * w;
    z = (D + omega * L_slt') \ v;
end