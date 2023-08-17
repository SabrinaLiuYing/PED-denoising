function A = FormMatrix(u, alpha)
    % set up
    n = size(u, 1);
    m = sqrt(n);
    h = 1 / (m + 1);
    u = reshape( u, [m, m])';
    beta = 1e-6;
    AW = zeros(m, m);
    AE = zeros(m, m);
    AS = zeros(m, m);
    AN = zeros(m, m);
    AC = zeros(m, m);

    % calcualte AW,AE,AS,AN,AC
    function val=inBound(i, j)
        if (i <= m && i >= 1 && j <= m && j >= 1)
            val = u(i, j);
        else
            val = 0;
        end
    end
    
    function val = inSquare(uTemp,u_2)
        val = ((uTemp-u_2)/h)^2;
    end

    function val = singleTerm(sTemp,s_2)
        val = 1 / (2 * (sTemp + s_2+ beta));
    end
    
    function val = total(tTemp,t_2)
        val = (-alpha / (h ^ 2)) * (tTemp+t_2);
    end

    for i = 1:m
        for j = 1:m
            AW(i, j) = total(singleTerm(inSquare(inBound(i,j),inBound(i - 1, j)),inSquare(inBound(i, j),inBound(i, j - 1))), ...
                singleTerm(inSquare(inBound(i , j),inBound(i - 1, j )),inSquare(inBound(i-1,j+1),inBound(i-1,j))));
            AE(i, j) = total(singleTerm(inSquare(inBound(i+1,j),inBound(i,j)),inSquare(inBound(i+1,j),inBound(i+1,j-1))), ...
                singleTerm(inSquare(inBound(i+1,j),inBound(i,j)),inSquare(inBound(i,j+1),inBound(i,j))));
            AS(i ,j) = total(singleTerm(inSquare(inBound(i,j),inBound(i-1,j)),inSquare(inBound(i,j),inBound(i,j-1))), ...
                singleTerm(inSquare(inBound(i+1,j-1),inBound(i,j-1)),inSquare(inBound(i,j),inBound(i,j-1))));
            AN(i, j) = total(singleTerm(inSquare(inBound(i+1,j),inBound(i,j)),inSquare(inBound(i,j+1),inBound(i,j))), ...
                singleTerm(inSquare(inBound(i,j+1),inBound(i-1,j+1)),inSquare(inBound(i,j+1),inBound(i,j))));
            AC(i, j) = - (AW(i, j) + AE(i, j) + AS(i, j) + AN(i, j)) + 1;
        end
    end

    ACTemp = AC';
    ANTemp = AN';
    ASTemp = AS';
    AETemp = AE(1:m - 1,:)';
    AWTemp = AW(2:m,:)';

    ACFinal = ACTemp(:);
    ANFinal = ANTemp(:);
    ASFinal = ASTemp(:);
    AEFinal = AETemp(:);
    AWFinal = AWTemp(:);

    ANFinal = ANFinal(1:n - 1);
    ASFinal = ASFinal(2:n);

    A = diag(ACFinal) + diag(ANFinal, 1) + diag(ASFinal, -1) + diag(AEFinal, m) + diag(AWFinal, -m);


    for i=1:m-1
        A(i * m, i * m + 1) = 0;
        A(i * m + 1, i * m) = 0;
    end

    A = sparse(A);
end
