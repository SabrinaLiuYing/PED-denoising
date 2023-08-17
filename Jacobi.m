function [x, iter] = Jacobi(A, b, x_initial, maxiter, tol)
 
    x = x_initial;

    for iter=0:maxiter
        r = (b - A * x);

        if (mod(iter, 1000) == 0)
            if norm(r) <= tol * norm(b)
                break
            end
        end

        x = x + diag(1./diag(A)) * r;
    end
end
