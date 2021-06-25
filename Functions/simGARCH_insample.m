
function out=simGARCH(rets,size,parameters);
Garch = NaN(size,1); 
Garch(1,1) = parameters(1,1)/(1-parameters(2,1)-parameters(3,1));
for j =2:(size+1)
    Garch(j,1) =   parameters(1,1)+ parameters(2,1)*(rets(j-1,1).^2)+ parameters(3,1)*Garch(j-1,1);
end

out = Garch(:,1);