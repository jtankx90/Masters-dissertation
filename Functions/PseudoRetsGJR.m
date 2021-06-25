function out = PseudoRetsGJR(residuals,samplesize,parameters)
pseudorets = NaN(samplesize,1);
repgjr = NaN(samplesize,1);
repgjr(1,1) = parameters(1,1)/(1-parameters(2,1)-parameters(4,1)-(parameters(3,1)/2));
pseudorets(1,1) = sqrt(repgjr(1,1))*datasample(residuals(:,1),1);

 for j=2:(samplesize)
    if pseudorets(j-1,1)>=0
        repgjr(j,1)= parameters(1,1) + (parameters(2,1)*(pseudorets(j-1,1).^2)) + (parameters(4,1)*repgjr(j-1,1));
    else
        repgjr(j,1) = parameters(1,1)+ (parameters(2,1)+parameters(3,1))*(pseudorets(j-1,1).^2)+ parameters(4,1)*repgjr(j-1,1);
    end    
    pseudorets(j,1) = sqrt(repgjr(j,1))*datasample(residuals(:,1),1);
 end   
      


out = pseudorets;