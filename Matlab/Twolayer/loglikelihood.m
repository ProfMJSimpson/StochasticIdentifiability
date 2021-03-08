function [ll] = loglikelihood(T, l, x0, p, type)

switch type
    case  'exact'
        f=ExitTime_ExactDiscreteQ(l,x0,p);
    case 'approx'
        f = ExitTime_ApproxGamma(l,x0,p);
end

ll=sum(log(f(T)));

end

