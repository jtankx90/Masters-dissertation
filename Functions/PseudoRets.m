function out = PseudoRets(residuals,samplesize,parameters)
pseudorets = NaN(samplesize,1);
repgarch = NaN(samplesize,1);
repgarch(1,1) = parameters(1,1)/(1-parameters(2,1)-parameters(3,1));
pseudorets(1,1) = sqrt(repgarch(1,1))*datasample(residuals(:,1),1);

 for j=2:(samplesize)
    repgarch(j,1)= parameters(1,1) + (parameters(2,1)*(pseudorets(j-1,1).^2)) + (parameters(3,1)*repgarch(j-1,1));
    pseudorets(j,1) = sqrt(repgarch(j,1))*datasample(residuals(:,1),1);
 end   
      


out = pseudorets;