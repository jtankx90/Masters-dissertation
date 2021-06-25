clear all

oridatasp500 = readtable('dissertation data-finalv2.xlsx','Sheet','SP500');
ret = [table2array(oridatasp500(:,3))*100];
dt = datetime(table2array(oridatasp500(:,1)),'InputFormat','dd/MM/yyyy');

len = [length(ret)];


%Note: require manual change of backtesting window size to 250 and 1250
%samplesize=len-1250
samplesize= len-250;
Size = len-samplesize+1;
index = [0.05,0.01];



vol = NaN(Size,4);
dof = NaN(Size,2);
var= NaN(Size,4);
res = NaN(samplesize,Size);
resT = NaN(samplesize,Size);
gjrRES = NaN(samplesize,Size);
TgjrRES =NaN(samplesize,Size);



tic
parfor i=1:Size
        display(i)
       %% initalise
       %GARCH
       window = ret(i:(samplesize+i-1));

       %estimation
      [Parameters,LL]= tarch(window,1,0,1);
      [GJRParameters,LL1] =tarch(window,1,1,1);
      [Tparameters,LL2] = tarch(window,1,0,1,'STUDENTST');
      [gjrTparams,LL3]=tarch(window,1,1,1,'STUDENTST');
      dof(i,:) = [Tparameters(4,1),gjrTparams(5,1)];


      %recursion of variance using function. A function is written to make
      %use of the parallel computing toolbox
      var(i,:) = [simGARCH(window,samplesize,Parameters),simGJR(window,samplesize,GJRParameters)...
                  ,simGARCH(window,samplesize,Tparameters),simGJR(window,samplesize,gjrTparams)];
      res(:,i) = simGARCHres(window,samplesize,Parameters);
      gjrRES(:,i) = simGJRres(window,samplesize,GJRParameters);
      resT(:,i) = simGARCHres(window,samplesize,Tparameters);
      TgjrRES(:,i) = simGJRres(window,samplesize,gjrTparams);
end 
toc

vol = sqrt(var(1:end-1,:));
retwindow = ret(samplesize+1:end);



%% VaR 
%student-t standard deviation
std                   = sqrt(minus(dof,2)./dof);


%5 percent VaR
vr5                   = NaN(1,4);
VaRest5               = NaN(Size-1,4);
VaRest5(:,1:2)        = vol(:,1:2).*-norminv(index(1)); %normal
VaRest5(:,3:4)        = vol(:,3:4).*std(1:end-1,:).*-tinv(index(1),dof(1:end-1,:)); %t-dist

%1-percent VaR
vr1                   =NaN(1,4);
VaRest1               =NaN(Size-1,4);
VaRest1(:,1:2)        =vol(:,1:2).*-norminv(index(2)); %normal
VaRest1(:,3:4)        =vol(:,3:4).*std(1:end-1,:).*-tinv(index(2),dof(1:end-1,:));  %t-dist


%% ES
% 5 percent ES
ES5(:,1:2)          =vol(:,1:2).*(normpdf(norminv(index(1)))/index(1)); %normal

%t-dist
A5                	=(tpdf(tinv(index(1),dof(1:end-1,:)),dof(1:end-1,:)))./(index(1)); 
B5                  =(dof(1:end-1,:)+ tinv(index(1),dof(1:end-1,:)).^2)./minus(dof(1:end-1,:),1);
ES5(:,3:4)          = vol(:,3:4).*std(1:end-1,:).*A5.*B5;



ES1(:,1:2)         =vol(:,1:2).*(normpdf(norminv(index(2)))/index(2)); %normal

%t-dist
A1                 =(tpdf(tinv(index(2),dof(1:end-1,:)),dof(1:end-1,:)))./(index(2));
B1                 =(dof(1:end-1,:)+ tinv(index(2),dof(1:end-1,:)).^2)./minus(dof(1:end-1,:),1);
ES1(:,3:4)         = vol(:,3:4).*std(1:end-1,:).*A1.*B1;




%% Backtesting VaR
%VR

Breach5 = NaN(length(vol),4);
Breach1=NaN(length(vol),4);
for i=1:4
    Breach5(:,i) = ret(samplesize+1:end) < -VaRest5(:,i);
    Breach1(:,i) = ret(samplesize+1:end) < -VaRest1(:,i);
end
vio5 = sum(Breach5);
vio1 = sum(Breach1);
vr5(:,:)= vio5./(index(1)*(Size-1));
vr1(:,:)= vio1./(index(2)*(Size-1));



% Unconditional Coverage Test
%5 percent
phat5(:,:) = vio5/length(Breach5);
V1_5 = NaN(1,4);
V0_5 = NaN(1,4);
V1_5(:,:) = vio5;
V0_5(:,:) = length(Breach5) - V1_5(:,:);

LogLConst(1,:) = V1_5(:,:).*log(index(1)) + V0_5(:,:).*log(1-index(1));
LogLUnconst(1,:) = V1_5(:,:).*log(phat5(:,:)) + V0_5(:,:).*log(1-phat5(:,:));
TestStatisticUC(1,:) = -2.*minus(LogLConst(1,:),LogLUnconst(1,:));
SignificanceUC(1,:) = minus(1,chi2cdf(TestStatisticUC(1,:), 1));


% 1 percent
phat1(:,:) = vio1/length(Breach1);
V1_1 = NaN(1,4);
V0_1 = NaN(1,4);
V1_1(:,:) = vio1;
V0_1(:,:) = length(Breach1)-V1_1(:,:);

LogLConst(2,:) = V1_1(:,:).*log(index(2)) + V0_1(:,:).*log(1-index(2));
LogLUnconst(2,:) = V1_1(:,:).*log(phat1(:,:)) + V0_1(:,:).*log(1-phat1(:,:));
TestStatisticUC(2,:) = -2.*minus(LogLConst(2,:),LogLUnconst(2,:));
SignificanceUC(2,:) = minus(1,chi2cdf(TestStatisticUC(2,:), 1));

% Independence Test
ind = NaN(2,4);
testsigcc = NaN(2,4);
for i =1:4;
    ind(1,i) = ind_test(Breach5(:,i));
    ind(2,i) = ind_test(Breach1(:,i));
    testsigcc(1,i) = 1-chi2cdf(ind(1,i),1);
    testsigcc(2,i) = 1-chi2cdf(ind(2,i),1);
end;

LRCC = NaN(2,4);
LRpvalues = NaN(2,4);

%cc
for i =1:4;
   LRCC(1,i) =   TestStatisticUC(1,i) + ind(1,i);
   LRCC(2,i) = TestStatisticUC(2,i) + ind(2,i);
   LRpvalues(1,i) = 1-chi2cdf(LRCC(1,i),2);
   LRpvalues(2,i) = 1-chi2cdf(LRCC(2,i),2);
end


% %loss function (Not in USE)
% %lossfun
% lossfunc5 = zeros(Size-1,4);
% for i =1:4
%     p = find(retwindow < -VaRest5(:,i));
%     lossfunc5(p,i) = (retwindow(p)- -VaRest5(p,i)).^2; 
% end
% loss5 = sum(lossfunc5)
% 
% 
% lossfunc1 = zeros(Size-1,4);
% for i =1:4
%     p = find(retwindow < -VaRest1(:,i));
%     lossfunc1(p,i) = (retwindow(p)- -VaRest1(p,i)).^2; 
% end
% loss1 = sum(lossfunc1)


%% Backtesting ES
combinedvar= [VaRest5, VaRest1];
combinedes = [ES5, ES1];
nES = NaN(2,4);

for i =1:4
    q = find(retwindow < -combinedvar(:,i));
    nES(1,i) = mean(retwindow(q) ./ -combinedes(q,i));
    q1 = find(retwindow <= -combinedvar(:,i+4));
    nES(2,i)= mean(retwindow(q1) ./ -combinedes(q1,i+4));
end





%% Plot
Figno = 1
models = {'GARCH(1,1)','GJR(1,1)','GARCH-T(1,1)','GJR-T(1,1)'}
% 5 percent
for i= 1:4
    subplot(2,2,[i])
    figure(Figno);
    chartdata = [retwindow,-vol(:,i),-VaRest5(:,i),-ES5(:,i)];
    TitleStr = strcat({'S&P Returns with '},models(i),{' volatility and 5% VaR '});
    title(TitleStr) 
    grid on;
    hold on;
    plot(dt(len-Size+2:end),chartdata(:,1), 'LineStyle', '--' , ...
                      'LineWidth', 0.5, 'Color', 'black');
    plot(dt(len-Size+2:end),chartdata(:,2), 'LineStyle', '-.', ...
                     'LineWidth', 1, 'Color', 'magenta');
    plot(dt(len-Size+2:end),chartdata(:,3), 'LineStyle', '-', ...
                     'LineWidth', 1, 'Color', 'red');   
    plot(dt(len-Size+2:end),chartdata(:,4), 'LineStyle', '-', ...
                     'LineWidth', 1, 'Color', 'blue');   
    legend(['S&P Returns',models(i),'VaR','ES']);
    legend('boxoff')
    set(gcf, 'color', 'white');
    
end


for i= 1:4
    subplot(2,2,[i])
    figure(Figno);
    chartdata = [retwindow,-vol(:,i),-VaRest1(:,i),-ES1(:,i)];
    TitleStr = strcat({'S&P Returns with '},models(i),{' volatility and 1% VaR '});
    title(TitleStr) 
    grid on;
    hold on;
    plot(dt(len-Size+2:end),chartdata(:,1), 'LineStyle', '--' , ...
                      'LineWidth', 0.5, 'Color', 'black');
    plot(dt(len-Size+2:end),chartdata(:,2), 'LineStyle', '-.', ...
                     'LineWidth', 1, 'Color', 'magenta');
    plot(dt(len-Size+2:end),chartdata(:,3), 'LineStyle', '-', ...
                     'LineWidth', 1, 'Color', 'red');   
    plot(dt(len-Size+2:end),chartdata(:,4), 'LineStyle', '-', ...
                     'LineWidth', 1, 'Color', 'blue');   
    legend(['S&P Returns',models(i),'VaR','ES']);
    legend('boxoff')
    set(gcf, 'color', 'white');
    
end



