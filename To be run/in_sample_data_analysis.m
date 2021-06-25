clear all
oridatasp500 = readtable('dissertation data-finalv2.xlsx','Sheet','SP500');
ret = [table2array(oridatasp500(:,3))*100];
len = length(ret);



%% Descriptive Statistics

%SP500 price and return plot
dt = datetime(table2array(oridatasp500(:,1)),'InputFormat','dd/MM/yyyy');
figure(3);
subplot(2,1,1)
TitleStr2 = strcat({'SP500 Prices'});
title(TitleStr2) 
hold on
price = [table2array(oridatasp500(:,2))];

plot(dt,price, 'LineStyle', '--' , ...
                  'LineWidth', 0.5, 'Color', 'blue');


subplot(2,1,2)
TitleStr = strcat({'SP500 Returns'});
title(TitleStr) 
hold on
plot(dt,ret/100, 'LineStyle', '--' , ...
                  'LineWidth', 0.5, 'Color', 'blue');
stockstats = NaN(6,1);

%divided by 100 to negate the multiplication earlier
stockstats(1,1) = mean(ret)/100;
stockstats(2,1) = std(ret)/100;
stockstats(3,1) = skewness(ret);
stockstats(4,1) = kurtosis(ret);

[~,P,JBSTAT] = jbtest(ret);
stockstats(5,1) = JBSTAT;
stockstats(6,1) = P;
% Ljung-Box Qtest
[h,p,stat]=lbqtest(ret.^2,'Lags',[1:8]);


% QQ Plot
subplot(2,1,1)
qqplot(ret)
subplot(2,1,2)
qqplot(ret, fitdist(ret,'tLocationScale'))

% Density Plot
y = normrnd(mean(ret),std(ret),[len,1]);
[t1,t2] = ksdensity(ret,'Function','pdf');
plot(t2,t1,'Linewidth',1);
hold on
[t3,t4]=ksdensity(y,'Function','pdf');
plot(t4,t3,'Linewidth',1);
%xlim([min(xlim),0.1]);
ylim([0,65]);
lgd = legend(['SP 500'],['Simulated Normal']);
lgd.FontSize = 12.5;
title('Empirical PDF plot of SP500 returns vs Random Normal');
hold off


% ACF 
figno=1;
figure(figno);
autocorr(ret);
hold off
Titlestr = ['Sample ACF for SP500 returns']
title(Titlestr);
ylabel('Historical Sample')
figno = figno+1;


figure(figno);
autocorr(ret.^2);
hold off
Titlestr = ['Sample ACF for SP500 squared returns']
title(Titlestr);
ylabel('Historical Sample')
figno = figno+1;


%% In-sample GARCH
len = length(ret);
vol = NaN(len,1);
var= NaN(len,4);
res = NaN(len,4);

%GARCH

%estimation
[Parameters,LL,~,VCVROBUST,~,score]= tarch(ret,1,0,1);
[GJRParameters,LL1,~,VCVROBUST1] =tarch(ret,1,1,1);
[Tparameters,LL2,~,VCVROBUST2] = tarch(ret,1,0,1,'STUDENTST');

[gjrTparams,LL3,~,VCVROBUST3]=tarch(ret,1,1,1,'STUDENTST');

var(:,:) = [simGARCH_insample(ret,len,Parameters),simGJR_insample(ret,len,GJRParameters)...
              ,simGARCH_insample(ret,len,Tparameters),simGJR_insample(ret,len,gjrTparams)];
res(:,:) = [simGARCHres(ret,len,Parameters),simGJRres(ret,len,GJRParameters),simGARCHres(ret,len,Tparameters)...
            ,simGJRres(ret,len,gjrTparams)];
vol = sqrt(var(1:end-1,1:end));
        
%coefficient test
normSTD = sqrt(diag(VCVROBUST));
gjrSTD = sqrt(diag(VCVROBUST1));
Tstd = sqrt(diag(VCVROBUST2));
gjrTstd = sqrt(diag(VCVROBUST3));

tstats = NaN(4,5);
tstats(1,1:length(normSTD)) = (Parameters./normSTD).'
tstats(2,1:length(gjrSTD))= (GJRParameters./gjrSTD).'
tstats(3,1:length(Tstd))= (Tparameters./Tstd).'
tstats(4,1:length(gjrTstd))= (gjrTparams./gjrTstd).'
tpvalue = 1-normcdf(tstats);

%Likelihood ratio test
logLL = [LL,LL1,LL2,LL3];
LR = 2*minus(logLL(:,2:4),logLL(:,1));
pvalues = 1-chi2cdf(LR,[1,1,2]);        
        






%% QQ plot residuals


Figure(1)
for i=1:2
    subplot(4,2,[i])  
    if i ==1
        qqplot(res(:,1));
        Title = strcat('QQ plot of   ', models(1), 'residuals against Normal')
        title(Title)
    else
        qqplot(res(:,1),fitdist(res(:,i),'tLocationScale'));
        Title = strcat('QQ plot of   ', models(1), 'residuals against Student-t')
        title(Title)
    end
end

hold on

for i=1:2
    subplot(4,2,[i]+2)  
    if i ==1
        qqplot(res(:,2));
        Title = strcat('QQ plot of   ', models(2), 'residuals against Normal')
        title(Title)
    else
        qqplot(res(:,2),fitdist(res(:,i),'tLocationScale'));
        Title = strcat('QQ plot of   ', models(2), 'residuals against Student-t')
        title(Title)
    end
end
hold on

for i=1:2
    subplot(4,2,[i]+4)  
    if i ==1
        qqplot(res(:,3));
        Title = strcat('QQ plot of   ', models(3), 'residuals against Normal')
        title(Title)
    else
        qqplot(res(:,3),fitdist(res(:,i),'tLocationScale'));
        Title = strcat('QQ plot of   ', models(3), 'residuals against Student-t')
        title(Title)
    end
end
hold on


for i=1:2
    subplot(4,2,[i]+6)  
    if i ==1
        qqplot(res(:,4));
        Title = strcat('QQ plot of   ', models(4), 'residuals against Normal')
        title(Title)
    else
        qqplot(res(:,4),fitdist(res(:,i),'tLocationScale'));
        Title = strcat('QQ plot of   ', models(4), 'residuals against Student-t')
        title(Title)
    end
end


%% ACF Plot of residuals

figure(2)
for i=1:2
    subplot(4,2,[i])  
    if i ==1
        autocorr(res(:,1));
        Title = strcat('ACF plot of   ', models(1), 'residuals')
        title(Title)
    else
        autocorr(res(:,1).^2);
        Title = strcat('ACF plot of  squared  ', models(1), 'residuals ')
        title(Title)
    end
end

hold on
for i=1:2
    subplot(4,2,[i]+2)  
    if i ==1
        autocorr(res(:,2));
        Title = strcat('ACF plot of   ', models(2), 'residuals')
        title(Title)
    else
        autocorr(res(:,2).^2);
        Title = strcat('ACF plot of  squared  ', models(2), 'residuals ')
        title(Title)
    end
end
hold on

for i=1:2
    subplot(4,2,[i]+4)  
    if i ==1
        autocorr(res(:,3));
        Title = strcat('ACF plot of   ', models(3), 'residuals')
        title(Title)
    else
        autocorr(res(:,3).^2);
        Title = strcat('ACF plot of  squared  ', models(3),  'residuals ')
        title(Title)
    end
end
hold on


for i=1:2
    subplot(4,2,[i]+6)  
    if i ==1
        autocorr(res(:,4));
        Title = strcat('ACF plot of  ', models(4), 'residuals')
        title(Title)
    else
        autocorr(res(:,4).^2);
        Title = strcat('ACF plot of  squared  ', models(4),  'residuals')
        title(Title)
    end
end


%residual normality and autocorrelation test
jbpvalues = NaN(1,4);
lbq = NaN(1,4);

for i=1:4
    [~,P2,JBSTAT]=jbtest(res(:,i));
    [h,p2,stat] = lbqtest(res(:,i).^2,'Lags',[8]);
    jbpvalues(1,i) = P2;
    lbq(1,i) = p2;
end;
 




