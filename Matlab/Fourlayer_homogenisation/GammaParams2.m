function [a,b] = GammaParams2(M,V)

    

    a = @(x0,p1,p2) M(x0,p1,p2)^2 / V(x0,p1,p2);
    b = @(x0,p1,p2) V(x0,p1,p2) / M(x0,p1,p2);

end