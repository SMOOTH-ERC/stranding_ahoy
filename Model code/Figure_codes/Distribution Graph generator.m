x = -2:0.01:2;
pd1 = makedist('Stable','alpha',1.1,'beta',-1,'gam',(sum(((beta.^(1:S).*sigma_u_H(t,t+1:t+S))*gamma_H(t)).^1.1))^(1/1.1),'delta',0)
pd2 = makedist('Stable','alpha',1.2,'beta',-1,'gam',(sum(((beta.^(1:S).*sigma_u_H(t,t+1:t+S))*gamma_H(t)).^1.2))^(1/1.2),'delta',0)
pd25 = makedist('Stable','alpha',1.5,'beta',-1,'gam',(sum(((beta.^(1:S).*sigma_u_H(t,t+1:t+S))*gamma_H(t)).^1.5))^(1/1.5),'delta',0)
pd3 = makedist('Stable','alpha',1.7,'beta',-1,'gam',(sum(((beta.^(1:S).*sigma_u_H(t,t+1:t+S))*gamma_H(t)).^1.7))^(1/1.7),'delta',0)
pd4 = makedist('Stable','alpha',1.9,'beta',-1,'gam',(sum(((beta.^(1:S).*sigma_u_H(t,t+1:t+S))*gamma_H(t)).^1.9))^(1/1.9),'delta',0)
pdf1 = pdf(pd1,x);
pdf2 = pdf(pd2,x);
pdf25 = pdf(pd25,x);
pdf3 = pdf(pd3,x);
pdf4 = pdf(pd4,x);
pdfnorm = normpdf(x, 0, sigma_r_H(t));
figure
plot(x,pdf1,'b--');
hold on
plot(x,pdf2,'r--');
plot(x,pdf25,'--');
plot(x,pdf3,'p--');
plot(x,pdf4,'g--');
plot(x, pdfnorm, 'k-')
title('Compare Alpha Parameters in Stable Distribution PDF Plots')
legend('\alpha = 1.1','\alpha = 1.5','\alpha = 1.7','\alpha = 1.9', 'Normal', 'Location','northwest')
hold off

max(finverse(pdf1))
Beta = beta*ones(S,1)
Beta
