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
var= NaN(Size,4);

%garch
resQuantilecentered5 = NaN(Size,round(samplesize*index(1)));
resQuantilecentered1 = NaN(Size,round(samplesize*index(2)));
%garch-t
resTQuantilecentered5 = NaN(Size,round(samplesize*index(1)));
resTQuantilecentered1 = NaN(Size,round(samplesize*index(2)));
%gjr
GJRresQuantilecentered5= NaN(Size,round(samplesize*index(1)));
GJRresQuantilecentered1 = NaN(Size,round(samplesize*index(2)));

%gjr-t
gjrTresQuantilecentered5 = NaN(Size,round(samplesize*index(1)));
gjrTresQuantilecentered1 = NaN(Size,round(samplesize*index(2)));

tic
parfor i= 1:(Size)
        display(i)
       %% initalise
       %GARCH estimation
       window = ret(i:samplesize+i-1);

      [Parameters]= tarch(window,1,0,1);
      [GJRParameters] =tarch(window,1,1,1);
      [Tparameters] = tarch(window,1,0,1,'STUDENTST');
      [gjrTparams]=tarch(window,1,1,1,'STUDENTST');




       % 1) retrieve residuals
       % 2) variance is initialise with unconditional variance computed from estimated parameters
       % 3) de-mean residuals & sort
      res = simGARCHres(window,samplesize,Parameters);
      centeredres = minus(res,mean(res));
      centeredres = sort(centeredres);

      gjrRES = simGJRres(window,samplesize,GJRParameters);
      centeredgjrRES = minus(gjrRES, mean(gjrRES));
      centeredgjrRES = sort(centeredgjrRES);


      resT = simGARCHres(window,samplesize,Tparameters);
      centeredresT = minus(resT, mean(resT));
      centeredresT = sort(centeredresT);


      TgjrRES = simGJRres(window,samplesize,gjrTparams);
      centeredTgjrRES = minus(TgjrRES,mean(TgjrRES));
      centeredTgjrRES = sort(centeredTgjrRES);



       %retrieve residuals for VaR and ES calculation



      %garch
      resQuantilecentered5(i,:) = centeredres((1:round(samplesize*index(1))),:);
      resQuantilecentered1(i,:) = centeredres((1:round(samplesize*index(2))),:);

      %gjr
      GJRresQuantilecentered5(i,:) = centeredgjrRES((1:round(samplesize*index(1))),:);
      GJRresQuantilecentered1(i,:) = centeredgjrRES((1:round(samplesize*index(2))),:)

      %garch-t
      resTQuantilecentered5(i,:) = centeredresT((1:round(samplesize*index(1))),:);
      resTQuantilecentered1(i,:) =centeredresT((1:round(samplesize*index(2))),:);

      %gjr-t
      gjrTresQuantilecentered5(i,:) = centeredTgjrRES((1:round(samplesize*index(1))),:);
      gjrTresQuantilecentered1(i,:) = centeredTgjrRES((1:round(samplesize*index(2))),:);



      var(i,:) = [simGARCH(window,samplesize,Parameters),simGJR(window,samplesize,GJRParameters)...
                  ,simGARCH(window,samplesize,Tparameters),simGJR(window,samplesize,gjrTparams)];
end    
toc

vol =sqrt(var(1:end-1,:));
retwindow = ret(samplesize+1:end);

%% VaR 
%5percentVaR
%VaR = average of all bootstrap VaRs where each column represent each moving window, each row represent bootstraps
% Transpose of the VaR matrix such as the rows are the moving window
vr5                   = NaN(1,4);
VaRest5               = NaN(Size-1,4);
combinedres5          =[resQuantilecentered5(1:end-1,end),GJRresQuantilecentered5(1:end-1,end),resTQuantilecentered5(1:end-1,end)...
                      ,gjrTresQuantilecentered5(1:end-1,end)];
VaRest5(:,:)          = vol.*-combinedres5;

%1-percent VaR
vr1                   =NaN(1,4);
VaRest1               =NaN(Size-1,4);
combinedres1          =[resQuantilecentered1(1:end-1,end),GJRresQuantilecentered1(1:end-1,end),resTQuantilecentered1(1:end-1,end)...
                      ,gjrTresQuantilecentered1(1:end-1,end)];
VaRest1(:,:)          =vol.*-combinedres1;


%% ES
% ES is computed as the average of residuals
EScombinedres5      = [mean(resQuantilecentered5(1:end-1,(1:end)),2), mean(GJRresQuantilecentered5(1:end-1,(1:end)),2), mean(resTQuantilecentered5(1:end-1,(1:end)),2)...
                    ,mean(gjrTresQuantilecentered5(1:end-1,(1:end)),2)];
EScombinedres1      = [mean(resQuantilecentered1(1:end-1,(1:end)),2), mean(GJRresQuantilecentered1(1:end-1,(1:end)),2), mean(resTQuantilecentered1(1:end-1,(1:end)),2)...
                    , mean(gjrTresQuantilecentered1(1:end-1,(1:end)),2)];
ES5                 =vol.*-EScombinedres5;
ES1                 =vol.*-EScombinedres1;


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
SignificanceUC(2,:) = minus(1,chi2cdf(TestStatisticUC(2,:), 1));% Independence Test

% Independence Test
ind = NaN(2,4);
testsigcc = NaN(2,4);
for i =1:4;
    ind(1,i) = ind_test(Breach5(:,i));
    ind(2,i) = ind_test(Breach1(:,i));
    testsigcc(1,i) = 1-chi2cdf(ind(1,i),1);
    testsigcc(2,i) = 1-chi2cdf(ind(2,i),1);
end;


%Conditional Coverage
LRCC = NaN(2,4);
LRpvalues = NaN(2,4);

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

%% Backtesting Expected Shortfall

combinedvar= [VaRest5, VaRest1];
combinedes = [ES5, ES1];
nES = NaN(2,4);

for i =1:4
    q = find(retwindow <= -combinedvar(:,i));
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
