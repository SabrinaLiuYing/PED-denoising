beta = 1e-6;
K = 10;

%For SOR, the optimal alphas I found at each resolution were:
%16: 1.87
%32: 1.94
%64: 1.96
%in the case where the initial guess is warm-started with the current
%solution, uk.

for grid_size = 3:3
    omega = 0;
    n = 16*2^(grid_size-1);
    fprintf('\n\nCurrent dimensions: %dx%d\n',n,n);
    if(grid_size==1)
        alpha = 4e-2;
        omega = 1.87;
    elseif(grid_size==2)
        alpha = 3e-2;
        omega = 1.94;
    elseif(grid_size==3)
        alpha = 1.5e-2;
        omega = 1.96; 
    else
        alpha = 1.2e-2;
        omega = 1.98 %I didn't test this one.
    end
    [u_exact, z] = set_image(n);
%     imshow(u_exact)
%     imshow(z)
    u0 = FormRHS(z);
    
    for method = 4:4
        soln = u0; %reset the solution
        fprintf('\n');
        switch(method)
            case 1
                fprintf('Jacobi:\n');
            case 2
                fprintf('Gauss-Seidel:\n');
            case 3
                fprintf('SOR:\n');
            case 4
                fprintf('ConjugateGradient:\n');
        end
        start_t = cputime;
        iteration_total = 0;
        for k = 0:K
            A = FormMatrix(soln, alpha, beta, n);
            iter = 0;
            switch(method)
                case 1
                    [soln,iter] = Jacobi(A, u0, soln, 100000,1e-2);
                case 2
                    [soln,iter] = GaussSeidel(A, u0, soln, 100000,1e-2);
                case 3
                    [soln,iter] = SOR(omega, A, u0, soln, 100000,1e-2);
                case 4
                    [soln,iter] = ConjugateGradient(A, u0, soln, 100000,1e-2);
            end
            fprintf('Step %d took %d iterations.\n', k, iter);
            iteration_total = iteration_total + iter;
        end
        end_t = cputime;
        fprintf('Resolution %d for method %d took %fs and  %d iters\n', n, method, end_t-start_t, iteration_total);
            
        result = x1to2(soln,n,n);
        figure(method);
        colormap(gray);
        imagesc(result);
    end
end

