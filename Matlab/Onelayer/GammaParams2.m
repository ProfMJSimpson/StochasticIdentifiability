function [a,b] = GammaParams2(M,V)

    

    a = @(x0,p1) M(x0,p1)^2 / V(x0,p1);
    b = @(x0,p1) V(x0,p1) / M(x0,p1);

end