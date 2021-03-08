function [M,V] = MomentsCentral(M1,M2)

    M = M1;
    V = @(x0,p1,p2,p3) M2(x0,p1,p2,p3) - M1(x0,p1,p2,p3)^2;

end