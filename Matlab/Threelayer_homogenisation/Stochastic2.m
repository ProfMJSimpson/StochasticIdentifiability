%Stochastic model two layers Mat Simpson September 2020
function [T] =Stochastic2(l1, l2, x0, p1, p2, MC, varargin)
if nargin==7
rng(varargin{1});
end

%l1 length of first layer
%l2 length of second layer
%x0 release point
%p1 diffusivity in first layer
%p2 diffusivity in second layer
%MC number of realisations

T=zeros(1,MC);
I1=l1; % Location of the first interface
I2=l1+l2; %Location of the second interface
LL=l1+l2; %Length

        for i = 1:MC
        t=0;
        tau=1;
       

        x = x0;

            if (x==0)
            t=0;
            T(i)=t;
            %T2(i)=t^2;
            end

            while x > 0
    
   
                if(x > 0 && x < I1)
                    R=rand;
                    if R < p1
                    x = x - 1;
                    elseif(R > p1 && R < (p1+p1))
                    x = x + 1;
                    elseif(R > (p1+p1) && R < 1)
                    x = x;
                    end
                t = t+tau;    
                end
                
                 if (x==I1)
                    R=rand;
                    if R < p1
                    x = x - 1;
                    elseif(R > p1 && R < (p1+p2))
                    x = x + 1;
                    elseif(R > (p1+p2) && R < 1)
                    x = x;
                    end
                       t = t+tau;
                 end
               
                if(x > I1 && x < I2)
                    R=rand;
                    if R < p2
                    x = x - 1;
                    elseif(R > p2 && R < (p2+p2))
                    x = x + 1;
                    elseif(R > (p2+p2) && R < 1)
                    x = x;
                    end
                t = t+tau;    
                end                
                
                
         
                
               if (x==LL)
                   R=rand;
                    if R < p2
                    x = x - 1;
                    elseif(R > p2 && R < 1)
                        x = x;
                    end
                t = t+tau;
               end
        T(i)=t;
        %T2(i)=t^2;
       
            end

        end
  
        
   
 
