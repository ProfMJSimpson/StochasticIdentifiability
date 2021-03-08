clearvars
p1=0.2;
p2=0.4;
P0=[0.1, 0.1];
lb=[0.05, 0.5];
ub=[0.05, 0.5];
MC=500;
x01=30;
x02=70;
l1=30;
l2=40;
[T1] =Stochastic2(l1, l2, x01, p1, p2, MC,1);
[T2] =Stochastic2(l1, l2, x02, p1, p2, MC,1);
fune=@(p) -loglikelihood(T1, [l1, l2], x01, p, 'exact')...
          -loglikelihood(T2, [l1, l2], x02, p, 'exact');
      
%[p_mle_e,nLLmin] = fmincon(fune,P0,[],[],[],[],lb,ub);
%LLmaxe = -nLLmin;

funa=@(p) -loglikelihood(T1, [l1, l2], x01, p, 'approx')...
          -loglikelihood(T2, [l1, l2], x02, p, 'approx');

%[p_mle_a,nLLmin] = fmincon(funa,P0,[],[],[],[],lb,ub);
%[p_mle_a,nLLmin] = fminsearch(funa,P0); %,[],[],[],[],lb,ub);
%LLmaxa = -nLLmin;

h1=0.05:0.01:0.50;
h2=0.05:0.01:0.50;
for i=1:length(h1)
    for j=1:length(h2)
    yy(i,j) = -fune([h1(i),h2(j)]); 
    zz(i,j) = -funa([h1(i),h2(j)]);
    end
end

figure
contourf(h2,h1,yy)
colorbar
title('Exact Log-likelihood')
figure
contourf(h2,h1,zz)
colorbar
title('Approximate Log-likelihood')

figure
contourf(h2,h1,-yy)
colorbar
title('Exact Negative Log-likelihood')
figure
contourf(h2,h1,-zz)
colorbar
title('Approximate Negative Log-likelihood')

%Manually optimize out the nuisance params

for i=1:length(h1)
   ph1e(i)=max(yy(i,:)); %maximse out h2
   ph1a(i)=max(zz(i,:)); %maximse out h2
end

for i=1:length(h2)
   ph2e(i)=max(yy(:,i)); %maximise out h1
   ph2a(i)=max(zz(:,i)); %maximise out h1
end

figure
subplot(1,2,1)
plot(h1,ph1e-max(ph1e),'r')
hold on
plot(h1,ph1a-max(ph1a),'g')
xlim([0.01, 0.50])
ylim([-3, 0])
subplot(1,2,2)
plot(h2,ph2e-max(ph2e),'r')
hold on
plot(h2,ph2a-max(ph2a),'g')
xlim([0.01, 0.50])
ylim([-3, 0])

% profiles, not sure if h1/h2 order correct...
figure
hold on
plot(h2,exp(max(yy,[],1)-max(yy,[],'all')),'r')
plot(h2,exp(max(zz,[],1)-max(zz,[],'all')),'g')
xlim([0.05, 0.50])

figure
hold on
plot(h1,exp(max(yy,[],2)-max(yy,[],'all')),'r')
plot(h1,exp(max(zz,[],2)-max(zz,[],'all')),'g')
xlim([0.05, 0.50])
