function [x, iter] = SOR(omega, A, b, x_initial, maxiter, tol)
    x = x_initial;
    D = diag(diag(A) / omega);
    L = tril(A);

    for iter=0:maxiter
        r = (b - A * x);

        if (mod(iter, 10) == 0)
            if norm(r) <= tol * norm(b)
                break
            end
        end
        x = x + (D + L) \ r;
    end
end
