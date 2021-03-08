%Stochastic model generating exit time data for 1D random walk
%Mat Simpson September 2020
clear all
locations=[10, 20, 30, 40, 50, 60, 70, 80, 90, 100]; %Number of positions where exit time data is calculated
MCMax=10000; %Number of MC realizations 
xloc=zeros(1,numel(locations)); %vector holding the locations at which we compare the ctm and discrete results
First_Raw_Moment=zeros(1,numel(locations));
Second_Raw_Moment=zeros(1,numel(locations));
T=zeros(1,MCMax);
T2=zeros(1,MCMax);

l1=30;
l2=70;
LL=l1+l2;
I1=l1;
I2=LL;
p1=0.2;
p2=0.3267;
for L = 1:numel(locations)  
    xloc(L)=locations(L);
    
        for i = 1:MCMax
        t=0;
        tau=1;
       

        x = xloc(L);

            if (x==0)
            t=0;
            T(i)=t;
            T2(i)=t^2;
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
        T2(i)=t^2;
       
            end

        end
  
        
   First_Raw_Moment(L)=mean(T);
   Second_Raw_Moment(L)=mean(T2);
   
    
end

D1=p1;
D2=p2;
x1=0:0.1:I1;
T1=(0.1e1/D1*l1+0.1e1/D1*l2).*x1-0.1e1/D1*x1.^2/0.2e1;
V1=(0.2e1/0.3e1/D1/D2*l2^3+0.2e1/0.3e1/D1^2*l1^3+0.2e1/D1^2*l1^2*l2+0.2e1/D1^2*l1*l2^2).*x1+(-0.1e1/D1^2*l1/0.3e1-0.1e1/D1^2*l2/0.3e1)*x1.^3+0.1e1/D1^2*x1.^4/0.12e2;

x2=I1:0.1:I2;
T2=-0.1e1/D2*l1^2/0.2e1-l1/D2*l2+0.1e1/D1*l1^2/0.2e1+l1/D1*l2+(0.1e1/D2*l1+0.1e1/D2*l2).*x2-0.1e1/D2*x2.^2/0.2e1;
V2=l1^4/D2^2/0.12e2+l1^3/D2^2*l2/0.3e1-0.2e1/0.3e1*l1/D2^2*l2^3-l1^4/D2/D1/0.2e1-0.2e1*l1^3/D2/D1*l2-0.2e1*l1^2/D2/D1*l2^2+0.2e1/0.3e1*l1/D2/D1*l2^3+0.5e1/0.12e2*l1^4/D1^2+0.5e1/0.3e1*l1^3/D1^2*l2+0.2e1*l1^2/D1^2*l2^2+(-0.1e1/D2^2*l1^3/0.3e1-0.1e1/D2^2*l1^2*l2+0.2e1/0.3e1/D2^2*l2^3+0.1e1/D2/D1*l1^3+0.3e1/D2/D1*l1^2*l2+0.2e1/D2/D1*l1*l2^2).*x2+(0.1e1/D2^2*l1^2/0.2e1+0.1e1/D2^2*l1*l2-0.1e1/D2*l1^2/D1/0.2e1-0.1e1/D2*l1/D1*l2)*x2.^2+(-0.1e1/D2^2*l1/0.3e1-0.1e1/D2^2*l2/0.3e1)*x2.^3+0.1e1/D2^2*x2.^4/0.12e2;


figure
plot(xloc,First_Raw_Moment,'ro')
hold on
plot(x1,T1,'b')
plot(x2,T2,'b')
xlabel('Position')
ylabel('First moment')
axis([0 LL 0 1.1*max(First_Raw_Moment)])

line([I1,I1],[ 0 1.1*max(First_Raw_Moment)])
line([I2,I2],[ 0 1.1*max(First_Raw_Moment)])



figure
plot(xloc,Second_Raw_Moment,'go')
hold on
plot(x1,V1,'b')
plot(x2,V2,'b')
xlabel('Position')
ylabel('Second moment')
axis([0 LL 0 1.1*max(Second_Raw_Moment)])
hold on
line([I1,I1],[ 0 1.1*max(Second_Raw_Moment)])
line([I2,I2],[ 0 1.1*max(Second_Raw_Moment)])


