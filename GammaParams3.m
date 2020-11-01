function [a,b] = GammaParams3(M,V)

    

    a = @(x0,p1,p2,p3) M(x0,p1,p2,p3)^2 / V(x0,p1,p2,p3);
    b = @(x0,p1,p2,p3) V(x0,p1,p2,p3) / M(x0,p1,p2,p3);
    
end