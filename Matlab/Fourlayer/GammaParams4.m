function [a,b] = GammaParams4(M,V)

   
    a = @(x0,p1,p2,p3,p4) M(x0,p1,p2,p3,p4)^2 / V(x0,p1,p2,p3,p4);
    b = @(x0,p1,p2,p3,p4) V(x0,p1,p2,p3,p4) / M(x0,p1,p2,p3,p4);
    
end