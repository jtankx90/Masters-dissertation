
function out=simGJR(rets,size,parameters);
gjr = NaN(size,1); 
gjr(1,1) = parameters(1,1)/(1-parameters(2,1)-parameters(4,1)-(parameters(3,1)/2));


for j =2:(size+1)
    if rets(j-1,1) >= 0 
        gjr(j,1) =  parameters(1,1)+ parameters(2,1)*(rets(j-1,1).^2)+ parameters(4,1)*gjr(j-1,1);
    else
        gjr(j,1) = parameters(1,1) + (parameters(2,1)+parameters(3,1))*(rets(j-1,1).^2)+parameters(4,1)*gjr(j-1,1);
    end
end

out = gjr(:,1);