clearvars
p1=0.25;
p2=0.25;
P0=[0.5, 0.5];
lb=[0.01, 0.01];
ub=[0.5, 0.5];
MC=5000;
x01=70;
x02=70;
l1=30;
l2=40;
[T1] =Stochastic2(l1, l2, x01, p1, p2, MC,1);
[T2] =Stochastic2(l1, l2, x02, p1, p2, MC,1);
fune=@(p) -loglikelihood(T1, [l1, l2], x01, p, 'exact')...
          -loglikelihood(T2, [l1, l2], x02, p, 'exact');
funa=@(p) -loglikelihood(T1, [l1, l2], x01, p, 'approx')...
          -loglikelihood(T2, [l1, l2], x02, p, 'approx');

P1=0.01:0.001:0.50;
for i=1:length(P1)
   yy(i) = fune(P1(i)); 
   zz(i) = funa(P1(i));
end
yy=yy-min(yy);
zz=zz-min(zz);

plot(P1,-yy,'r')
hold on
plot(P1,-zz,'g')
xline(p1)

