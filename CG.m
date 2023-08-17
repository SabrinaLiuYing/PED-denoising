function [x, iter] = CG(A, b, x_initial, maxiter, tol)
    x = x_initial;
    g = A' * (A * x) - (A' * b);
    temp = g;

    for iter=1:maxiter
        if norm(g) <= tol
            break
        end

        sd = A'* (A * temp);
        g = g - (g' * temp)/(temp' * sd) * sd;
        x = x - (g' * temp)/(temp' * sd) * temp;
        beta = (g' * sd) / (temp' * sd);
        temp = g - beta * temp;
    end
end
