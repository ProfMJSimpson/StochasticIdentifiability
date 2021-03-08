function [M,V] = MomentsCentral(M1,M2)

    M = M1;
    V = @(x0,p1) M2(x0,p1) - M1(x0,p1)^2;

end