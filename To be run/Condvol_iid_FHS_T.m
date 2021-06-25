
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



B= 1000;
%boot variance array initialise
bootNormalGarchvar = NaN(B,Size);
bootTGarchvar= NaN(B,Size);
bootGJRvar = NaN(B,Size);
BootTGJRvar = NaN(B,Size);



%Residual for VaR
%garch
resQuantilecentered5 = NaN(B,Size);
resQuantilecentered1 = NaN(B,Size);
%garch-t
resTQuantilecentered5 = NaN(B,Size);
resTQuantilecentered1 = NaN(B,Size);
%gjr
GJRresQuantilecentered5= NaN(B,Size);
GJRresQuantilecentered1 = NaN(B,Size);
%gjr-t
gjrTresQuantilecentered5 = NaN(B,Size);
gjrTresQuantilecentered1 = NaN(B,Size);

%Residual for ES
%garch
resES5 = NaN(B,Size);
resES1 = NaN(B,Size);
%garch-t
resTES5= NaN(B,Size);
resTES1 = NaN(B,Size);
%gjr
GJRresES5 =NaN(B,Size);
GJRresES1 = NaN(B,Size);
%gjr-t
gjrTresES5 = NaN(B,Size);
gjrTresES1 = NaN(B,Size);
ret = ret(999:3018);
tic

parfor i= 1:(Size)
    %% Initial Estimatiion (Step 1)

    window = ret(i:samplesize+i-1);
    [Parameters]= tarch(window,1,0,1);
    [GJRParameters] =tarch(window,1,1,1);
    [Tparameters] = tarch(window,1,0,1,'STUDENTST');
    [gjrTparams]=tarch(window,1,1,1,'STUDENTST');

    %residuals
    res= simGARCHres(window,samplesize,Parameters);
    gjrRES = simGJRres(window,samplesize,GJRParameters);
    resT = simGARCHres(window,samplesize,Tparameters);
    TgjrRES = simGJRres(window,samplesize,gjrTparams);


    %de-mean residuals
    demeanres = minus(res,mean(res,1));
    demeangjr = minus(gjrRES,mean(gjrRES,1));
    demeanresT=  minus(resT,mean(resT,1));
    demeangjrT = minus(TgjrRES,mean(TgjrRES,1));


    %% Step 2
    % generate pseudo returns from i.i.d resampling centered residuals,
    % where B = number of bootstrap estimates

    for k=1:B
    % 1) Generate pseudo variance and pseudo returns from Parameters estimated earlier
    % 2) Initial pseudo variance = unconditional variance
    % 3) recursion of variance with pseudo returns 
    % 4) pseudo returns = pseudo volatility * i.i.d resample from centered residuals
    % 5) pseudo variance for GJR models: 
      %  if pseudo returns(j-1) >0
            % alpha * pseudo returns(j-1)
      % if pseudo returns(j-1) <0
            % alpha+gamma * pseudo returns(j-1)
        pseudorets = PseudoRets(demeanres,samplesize,Parameters);
        pseudoretsGJR = PseudoRetsGJR(demeangjr,samplesize,GJRParameters);
        pseudoretsT = PseudoRets(demeanresT,samplesize,Tparameters);
        pseudoretsGJRT = PseudoRetsGJR(demeangjrT,samplesize,gjrTparams);

        %re-estimate garch using pseudo returns to obtain bootstrap parameters
       [hatParameters]= tarch(pseudorets(:,1),1,0,1);
       [hatGJR]= tarch(pseudoretsGJR(:,1),1,1,1);
       [hatTparams]= tarch(pseudoretsT(:,1),1,0,1,'STUDENTST');
       [hatGJRTparams]= tarch(pseudoretsGJRT(:,1),1,1,1,'STUDENTST');

        %re-trieve bootstrap residuals from pseudo returns, pseudo variance is initialise with the bootstrap parameters
        bootres =simGARCHres(pseudorets,samplesize,hatParameters);
        bootgjrRES = simGJRres(pseudoretsGJR,samplesize,hatGJR);       
        bootresT = simGARCHres(pseudoretsT,samplesize,hatTparams);
        bootTgjrRES = simGJRres(pseudoretsGJRT,samplesize,hatGJRTparams);


        %de-mean residuals
        bootres = minus(bootres,mean(bootres,1));
        bootres = sort(bootres);

        bootgjrRES = minus(bootgjrRES,mean(bootgjrRES,1));
        bootgjrRES = sort(bootgjrRES);

        bootresT= minus(bootresT,mean(bootresT,1));
        bootresT= sort(bootresT);

        bootTgjrRES= minus(bootTgjrRES,mean(bootTgjrRES,1));
        bootTgjrRES = sort(bootTgjrRES);



        %5 Percent res for VaR & ES

        % VaR residual is the cut-off point; ES residual is the average of residuals from 1:cut-off point
        resQuantilecentered5(k,i) = bootres(round(samplesize*index(1)),:);
        GJRresQuantilecentered5(k,i) = bootgjrRES(round(samplesize*index(1)),:);
        resTQuantilecentered5(k,i) = bootresT(round(samplesize*index(1)),:);
        gjrTresQuantilecentered5(k,i)= bootTgjrRES(round(samplesize*index(1)),:);

        resES5(k,i) = mean(bootres(1:(round(samplesize*index(1))),:),1);
        GJRresES5(k,i) = mean(bootgjrRES(1:(round(samplesize*index(1))),:),1);     
        resTES5(k,i)= mean(bootresT(1:(round(samplesize*index(1))),:),1);
        gjrTresES5(k,i) = mean(bootTgjrRES(1:(round(samplesize*index(1))),:),1);


        %1 Percent res for VaR & ES
        resQuantilecentered1(k,i) = bootres(round(samplesize*index(2)),:);
        GJRresQuantilecentered1(k,i) = bootgjrRES(round(samplesize*index(2)),:);
        resTQuantilecentered1(k,i) = bootresT(round(samplesize*index(2)),:); 
        gjrTresQuantilecentered1(k,i)= bootTgjrRES(round(samplesize*index(2)),:);

        resES1(k,i) = mean(bootres(1:(round(samplesize*index(2))),:),1);
        GJRresES1(k,i) = mean(bootgjrRES(1:(round(samplesize*index(2))),:),1);     
        resTES1(k,i)= mean(bootresT(1:(round(samplesize*index(2))),:),1);
        gjrTresES1(k,i) = mean(bootTgjrRES(1:(round(samplesize*index(2))),:),1);


        % 1) initial bootstrap variance = unconditional variance from bootstrap parameters
        % 2) recursion of variance with actual returns and bootstrap parameters

        bootNormalGarchvar(k,i)= simGARCH(window,samplesize,hatParameters);
        bootGJRvar(k,i) = simGJR(window,samplesize,hatGJR);
        bootTGarchvar(k,i) = simGARCH(window,samplesize,hatTparams);
        BootTGJRvar(k,i) = simGJR(window,samplesize,hatGJRTparams);
    end
end

toc


%Excluding the last forecast
normvol = sqrt(bootNormalGarchvar(:,1:end-1));
gjrvol = sqrt(bootGJRvar(:,1:end-1));
tvol = sqrt(bootTGarchvar(:,1:end-1));
gjrTvol = sqrt(BootTGJRvar(:,1:end-1));

PI = 0.05;
bound = [round(PI/2*B),round((1-PI/2)*B)];

%% VaR (element by element multiplication)


%5percentVaR
%VaR = average of all bootstrap VaRs where each column represent each moving window, each row represent bootstraps
% Transpose of the VaR matrix such as the rows are the moving window
BootVaR5(:,1:2) = [mean(normvol.*-resQuantilecentered5(:,1:end-1),1).',mean(gjrvol.*-GJRresQuantilecentered5(:,1:end-1),1).']; 
BootVaR5(:,3:4) = [mean(tvol.*-resTQuantilecentered5(:,1:end-1),1).',mean(gjrTvol.*-gjrTresQuantilecentered5(:,1:end-1),1).'];

%1percentVaR
BootVaR1(:,1:2) = [mean(normvol.*-resQuantilecentered1(:,1:end-1),1).',mean(gjrvol.*-GJRresQuantilecentered1(:,1:end-1),1).']; 
BootVaR1(:,3:4) = [mean(tvol.*-resTQuantilecentered1(:,1:end-1),1).',mean(gjrTvol.*-gjrTresQuantilecentered1(:,1:end-1),1).'];



%% ES(element by element multiplication)
% 5 percent ES
BootES5(:,1:2) =  [mean(normvol.*-resES5(:,1:end-1),1).',mean(gjrvol.*-GJRresES5(:,1:end-1),1).'];
BootES5(:,3:4) =  [mean(tvol.*-resTES5(:,1:end-1),1).',mean(gjrTvol.*-gjrTresES5(:,1:end-1),1).'];


%1 percent ES
BootES1(:,1:2) =  [mean(normvol.*-resES1(:,1:end-1),1).',mean(gjrvol.*-GJRresES1(:,1:end-1),1).'];
BootES1(:,3:4) =  [mean(tvol.*-resTES1(:,1:end-1),1).',mean(gjrTvol.*-gjrTresES1(:,1:end-1),1).'];




%% Prediction Interval
% For each moving window, there are B VaR estimates. Sort B VaR estimates
% and take the lower & upper bound as alpha/2 and 1-(alpha/2)
% Transpose such that rows are the moving window, columns are the lower and upper bound

% 5 percent VaR
normPI5 = sort(normvol.*-resQuantilecentered5(:,1:end-1),1);   
normBootPI5 = normPI5(bound,:).'; 

gjrPI5 = sort(gjrvol.*-GJRresQuantilecentered5(:,1:end-1),1);
gjrBootPI5 = gjrPI5(bound,:).';


tPI5 = sort(tvol.*-resTQuantilecentered5(:,1:end-1),1);
tBootPI5 = tPI5(bound,:).';

tgjrPI5 = sort(gjrTvol.*-gjrTresQuantilecentered5(:,1:end-1),1);
tgjrBootPI5 = tgjrPI5(bound,:).';

combinedVaRPI5 = [normBootPI5,gjrBootPI5,tBootPI5,tgjrBootPI5];


% 1 percent VaR
normPI1 = sort(normvol.*-resQuantilecentered1(:,1:end-1),1);   % Sort VaRs
normBootPI1 = normPI1(bound,:).'; 

gjrPI1 = sort(gjrvol.*-GJRresQuantilecentered1(:,1:end-1),1);
gjrBootPI1 = gjrPI1(bound,:).';


tPI1 = sort(tvol.*-resTQuantilecentered1(:,1:end-1),1);
tBootPI1 = tPI1(bound,:).';

tgjrPI1 = sort(gjrTvol.*-gjrTresQuantilecentered1(:,1:end-1),1);
tgjrBootPI1 = tgjrPI1(bound,:).';

combinedVaRPI1 = [normBootPI1,gjrBootPI1,tBootPI1,tgjrBootPI1];





% Expected Shortfall

% 5percent ES
normESPI5 = sort(normvol.*-resES5(:,1:end-1),1);
normBootESPI5 = normESPI5(bound,:).';

gjrESPI5 = sort(gjrvol.*-GJRresES5(:,1:end-1),1);
gjrBootESPI5 = gjrESPI5(bound,:).';

tESPI5 = sort(tvol.*-resTES5(:,1:end-1),1);
tBootESPI5 = tESPI5(bound,:).';

tgjrESPI5 = sort(gjrTvol.*-gjrTresES5(:,1:end-1),1);
tgjrBootESPI5 = tgjrESPI5(bound,:).';

combinedESPI5 = [normBootESPI5,gjrBootESPI5,tBootESPI5,tgjrBootESPI5];



%1 percent ES
normESPI1 = sort(normvol.*-resES1(:,1:end-1),1);
normBootESPI1 = normESPI1(bound,:).';

gjrESPI1 = sort(gjrvol.*-GJRresES1(:,1:end-1),1);
gjrBootESPI1 = gjrESPI1(bound,:).';

tESPI1 = sort(tvol.*-resTES1(:,1:end-1),1);
tBootESPI1 = tESPI1(bound,:).';

tgjrESPI1 = sort(gjrTvol.*-gjrTresES1(:,1:end-1),1);
tgjrBootESPI1 = tgjrESPI1(bound,:).';

combinedESPI1 = [normBootESPI1,gjrBootESPI1,tBootESPI1,tgjrBootESPI1];


%avg PI across estimation window
avgVarPI5 = mean(combinedVaRPI5,1);
avgVarPI1 = mean(combinedVaRPI1,1);
avgESPI5 = mean(combinedESPI5,1);
avgESPI1 = mean(combinedESPI1,1);


%Width of VaR interval
varPIlen5 = NaN(Size-1,4);
varPIlen1 = NaN(Size-1,4);

for i=1:4
    varPIlen5(:,i) = minus(combinedVaRPI5(:,i*2),combinedVaRPI5(:,i*2-1));
    varPIlen1(:,i) =minus(combinedVaRPI1(:,i*2),combinedVaRPI1(:,i*2-1));
end

%Width of ES interval
esPIlen5 = NaN(Size-1,4);
esPIlen1 = NaN(Size-1,4);

for i=1:4
    esPIlen5(:,i) = minus(combinedESPI5(:,i*2),combinedESPI5(:,i*2-1));
    esPIlen1(:,i) =minus(combinedESPI1(:,i*2),combinedESPI1(:,i*2-1));
end


%Standard Error

stdVaR5 = [std(normPI5,1).',std(gjrPI5,1).',std(tPI5,1).',std(tgjrPI5,1).'];
stdVaR1 = [std(normPI1,1).',std(gjrPI1,1).',std(tPI1,1).',std(tgjrPI1,1).'];
stdES5 = [std(normESPI5,1).',std(gjrESPI5,1).',std(tESPI5,1).',std(tgjrESPI5,1).'];
stdES1 = [std(normESPI1,1).',std(gjrESPI1,1).',std(tESPI1,1).',std(tgjrESPI1,1).'];


avgSTDvar5 = mean(stdVaR5);
avgSTDvar1 = mean(stdVaR1);
avgSTDes5 = mean(stdES5);
avgSTDes1 = mean(stdES1);


% %Tabulised results
% Headers = {'Lower Bound', 'Upper Bound','% Width', 'SE'}
% varPIresults5 = [avgVarPI5(1,[1,3,5,7]);avgVarPI5(1,[2,4,6,8]);avgpercenwidth(1,1:4);avgSTDvar5];
% esPIresults5 = [avgESPI5(1,[1,3,5,7]);avgESPI5(1,[2,4,6,8]);avgpercenwidth(2,1:4);avgSTDes5];
% varPIresults1 = [avgVarPI1(1,[1,3,5,7]);avgVarPI1(1,[2,4,6,8]);avgpercenwidth(1,5:8);avgSTDvar1];
% esPIresults1 = [avgESPI1(1,[1,3,5,7]);avgESPI1(1,[2,4,6,8]);avgpercenwidth(2,5:8);avgSTDes1];








%% Backtesting VaR
%VR
Breach5 = NaN(Size-1,4);
Breach1 =NaN(Size-1,4);
for i=1:4
    Breach5(:,i) = ret((2020-Size+2):end) < -BootVaR5(:,i);
    Breach1(:,i) = ret((2020-Size+2):end) < -BootVaR1(:,i);
end
vio5 = sum(Breach5);
vio1 = sum(Breach1);
vr5(:,:)= vio5./(index(1)*(Size-1));
vr1(:,:)= vio1./(index(2)*(Size-1));


% Unconditional Coverage Test
% 5 percent
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

% cc test
LRCC = NaN(2,4);
LRpvalues= NaN(2,4);

for i =1:4;
   LRCC(1,i) =   TestStatisticUC(1,i) + ind(1,i);
   LRpvalues(1,i) = 1-chi2cdf(LRCC(1,i),2);
   LRCC(2,i) =   TestStatisticUC(2,i) + ind(2,i);
   LRpvalues(2,i) = 1-chi2cdf(LRCC(2,i),2);
end




%% Backtesting Expected Shortfall

combinedvar = [BootVaR5, BootVaR1];
combinedes = [BootES5,BootES1];
retwindow = ret(samplesize+2:end);
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
combinedvol = [mean(normvol,1).',mean(gjrvol,1).',mean(tvol,1).',mean(gjrTvol,1).'];
% 5 percent
for i= 1:4
    subplot(2,2,[i])
    figure(Figno);
    chartdata = [retwindow,-combinedvol(:,i),-BootVaR5(:,i),-BootES5(:,i)];
    TitleStr = strcat({'S&P Returns with '},models(i),{' volatility and 5% VaR and ES '});
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
    chartdata = [retwindow,-combinedvol(:,i),-BootVaR1(:,i),-BootES1(:,i)];
    TitleStr = strcat({'S&P Returns with '},models(i),{' volatility and 1% VaR and ES '});
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

%% Width Difference plot
%(The width plot can be plotted using either condvol_iid_T or condvol_iid_FHS_T...
                % Have to run condvol_iid_T first to save the data, followed by loading the data for plot)

lenchartdata1fhs= [retwindow,minus(varPIlen1(:,:),esPIlen1(:,:))];
lenchartdata5fhs = [retwindow,minus(varPIlen5(:,:),esPIlen5(:,:))];

save('FHSlenchartdata1.mat','FHSlenchartdata1')
save('FHSlenchartdata5.mat','FHSlenchartdata5')

load('lenchartdata1.mat')
load('lenchartdata5.mat')

%1% plot
figure(1)
subplot(2,1,1)
TitleStr = 'SP500 returns with width difference between 1% iid-Bootstrap VaR and ES';
title(TitleStr) 
grid on;
hold on;
plot(dt(len-Size+2:end),iidlenchartdata1(:,1), 'LineStyle', '-' , ...
                      'LineWidth', 1, 'Color', 'cyan');
plot(dt(len-Size+2:end),iidlenchartdata1(:,2), 'LineStyle', '-.', ...
                     'LineWidth', 2, 'Color', 'red');
plot(dt(len-Size+2:end),iidlenchartdata1(:,3), 'LineStyle', '-.', ...
                     'LineWidth', 2, 'Color', 'blue');
plot(dt(len-Size+2:end),iidlenchartdata1(:,4), 'LineStyle', '-.', ...
                     'LineWidth', 2, 'Color', 'magenta');
plot(dt(len-Size+2:end),iidlenchartdata1(:,5), 'LineStyle', '-.', ...
                     'LineWidth', 2, 'Color', 'black');          
legend('Sp500 Ret','GARCH','GJR','T-GARCH','T-GJR');
legend('boxoff')
set(gcf, 'color', 'white');


%FHS plot
subplot(2,1,2)
%figure(2)
TitleStr = 'SP500 returns with width difference between 1% iid-FHS VaR and ES';
title(TitleStr) 
grid on;
hold on;
plot(dt(len-Size+2:end),FHSlenchartdata1(:,1), 'LineStyle', '-' , ...
                      'LineWidth', 1, 'Color', 'cyan');
plot(dt(len-Size+2:end),FHSlenchartdata1(:,2), 'LineStyle', '-.', ...
                     'LineWidth', 2, 'Color', 'red');
plot(dt(len-Size+2:end),FHSlenchartdata1(:,3), 'LineStyle', '-.', ...
                     'LineWidth', 2, 'Color', 'blue');
plot(dt(len-Size+2:end),FHSlenchartdata1(:,4), 'LineStyle', '-.', ...
                      'LineWidth', 2, 'Color', 'magenta');
plot(dt(len-Size+2:end),FHSlenchartdata1(:,5), 'LineStyle', '-.', ...
                      'LineWidth', 2, 'Color', 'black');                 
legend('Sp500 Ret','GARCH','GJR','T-GARCH','T-GJR');
legend('boxoff')
set(gcf, 'color', 'white');


% 5% plot

figure(2)
subplot(2,1,1)
TitleStr = 'SP500 returns with width difference between 5% iid VaR and ES';
title(TitleStr) 
grid on;
hold on;
plot(dt(3018-Size+2:end),iidlenchartdata5(:,1), 'LineStyle', '-' , ...
                      'LineWidth', 1, 'Color', 'cyan');
plot(dt(3018-Size+2:end),iidlenchartdata5(:,2), 'LineStyle', '-.', ...
                     'LineWidth', 2, 'Color', 'red');
plot(dt(3018-Size+2:end),iidlenchartdata5(:,3), 'LineStyle', '-.', ...
                     'LineWidth', 2, 'Color', 'blue');
plot(dt(3018-Size+2:end),iidlenchartdata5(:,4), 'LineStyle', '-.', ...
                     'LineWidth', 2, 'Color', 'magenta');
plot(dt(3018-Size+2:end),iidlenchartdata5(:,5), 'LineStyle', '-.', ...
                     'LineWidth', 2, 'Color', 'black');          
legend('Sp500 Ret','GARCH','GJR','T-GARCH','T-GJR');
legend('boxoff')
set(gcf, 'color', 'white');


%FHS plot
subplot(2,1,2)
TitleStr = 'SP500 returns with width difference between 5% iid-FHS VaR and ES';
title(TitleStr) 
grid on;
hold on;
plot(dt(3018-Size+2:end),FHSlenchartdata5(:,1), 'LineStyle', '-' , ...
                      'LineWidth', 1, 'Color', 'cyan');
plot(dt(3018-Size+2:end),FHSlenchartdata5(:,2), 'LineStyle', '-.', ...
                     'LineWidth',2, 'Color', 'red');
plot(dt(3018-Size+2:end),FHSlenchartdata5(:,3), 'LineStyle', '-.', ...
                     'LineWidth', 2, 'Color', 'blue');
plot(dt(3018-Size+2:end),FHSlenchartdata5(:,4), 'LineStyle', '-.', ...
                     'LineWidth', 2, 'Color', 'magenta');
plot(dt(3018-Size+2:end),FHSlenchartdata5(:,5), 'LineStyle', '-.', ...
                     'LineWidth', 2, 'Color', 'black');                 
legend('Sp500 Ret','GARCH','GJR','T-GARCH','T-GJR');
legend('boxoff')
set(gcf, 'color', 'white');




%% width plot 
%(The width difference plot can be plotted using either condvol_iid_T or condvol_iid_FHS_T...
%Have to run condvol_iid_FHS_T first to save the data, followed by loading the data for plot)

widthchartdata1fhs = [retwindow,varPIlen1(:,:),esPIlen1(:,:)];
widthchartdata5fhs = [retwindow,varPIlen5(:,:),esPIlen5(:,:)];


save('widthchartdata1fhs.mat','widthchartdata1fhs')
save('widthchartdata5fhs.mat','widthchartdata5fhs')

load('iidwidthchartdata1.mat');
load(';iidwidthchartdata5.mat');

%1%
%Normal
figure(1)
subplot(2,2,1)
TitleStr = 'SP500 returns with 1% iid-Bootstrap GARCH VaR and ES interval width';
title(TitleStr) 
grid on;
hold on;

plot(dt(len-Size+2:end),iidwidthchartdata1(:,1), 'LineStyle', '--' , ...
                      'LineWidth', 1, 'Color', 'cyan');
plot(dt(len-Size+2:end),iidwidthchartdata1(:,2), 'LineStyle', '-.', ...
                     'LineWidth', 2, 'Color', 'red');
plot(dt(len-Size+2:end),iidwidthchartdata1(:,6), 'LineStyle', '-.', ...
                     'LineWidth', 2, 'Color', 'blue');     
legend('Sp500 Ret','VaR Interval interval','ES Interval Width');
legend('boxoff')
set(gcf, 'color', 'white');


%FHS plot
subplot(2,2,2)
TitleStr = 'SP500 returns with 1% iid-FHS GARCH VaR and ES interval width';
title(TitleStr) 
grid on;
hold on;
plot(dt(len-Size+2:end),widthchartdata1fhs(:,1), 'LineStyle', '--' , ...
                      'LineWidth', 1, 'Color', 'cyan');
plot(dt(len-Size+2:end),widthchartdata1fhs(:,2), 'LineStyle', '-.', ...
                     'LineWidth', 2, 'Color', 'red');
plot(dt(len-Size+2:end),widthchartdata1fhs(:,6), 'LineStyle', '-.', ...
                     'LineWidth', 2, 'Color', 'blue');              
legend('Sp500 Ret','VaR Interval interval','ES Interval Width');
legend('boxoff')
set(gcf, 'color', 'white');

subplot(2,2,3)
TitleStr = 'SP500 returns with 1% iid-Bootstrap GJR VaR and ES interval width';
title(TitleStr) 
grid on;
hold on;

plot(dt(len-Size+2:end),iidwidthchartdata1(:,1), 'LineStyle', '--' , ...
                      'LineWidth', 1, 'Color', 'cyan');
plot(dt(len-Size+2:end),iidwidthchartdata1(:,3), 'LineStyle', '-.', ...
                     'LineWidth', 2, 'Color', 'red');
plot(dt(len-Size+2:end),iidwidthchartdata1(:,7), 'LineStyle', '-.', ...
                     'LineWidth', 2, 'Color', 'blue');     
legend('Sp500 Ret','VaR Interval interval','ES Interval Width');
legend('boxoff')
set(gcf, 'color', 'white');


%FHS plot
subplot(2,2,4)
TitleStr = 'SP500 returns with 1% iid-FHS GJR VaR and ES interval width';
title(TitleStr) 
grid on;
hold on;
plot(dt(len-Size+2:end),widthchartdata1fhs(:,1), 'LineStyle', '--' , ...
                      'LineWidth', 1, 'Color', 'cyan');
plot(dt(len-Size+2:end),widthchartdata1fhs(:,3), 'LineStyle', '-.', ...
                     'LineWidth', 2, 'Color', 'red');
plot(dt(len-Size+2:end),widthchartdata1fhs(:,7), 'LineStyle', '-.', ...
                     'LineWidth', 2, 'Color', 'blue');              
legend('Sp500 Ret','VaR Interval interval','ES Interval Width');
legend('boxoff')
set(gcf, 'color', 'white');



%student-t
subplot(2,2,1)
TitleStr = 'SP500 returns with 1% iid-Bootstrap tGARCH VaR and ES interval width';
title(TitleStr) 
grid on;
hold on;

plot(dt(len-Size+2:end),iidwidthchartdata1(:,1), 'LineStyle', '--' , ...
                      'LineWidth', 1, 'Color', 'cyan');
plot(dt(len-Size+2:end),iidwidthchartdata1(:,4), 'LineStyle', '-.', ...
                     'LineWidth', 2, 'Color', 'red');
plot(dt(len-Size+2:end),iidwidthchartdata1(:,8), 'LineStyle', '-.', ...
                     'LineWidth', 2, 'Color', 'blue');     
legend('Sp500 Ret','VaR Interval interval','ES Interval Width');
legend('boxoff')
set(gcf, 'color', 'white');


%FHS plot
subplot(2,2,2)
TitleStr = 'SP500 returns with 1% iid-FHS tGARCH VaR and ES interval width';
title(TitleStr) 
grid on;
hold on;
plot(dt(len-Size+2:end),widthchartdata1fhs(:,1), 'LineStyle', '--' , ...
                      'LineWidth', 1, 'Color', 'cyan');
plot(dt(len-Size+2:end),widthchartdata1fhs(:,4), 'LineStyle', '-.', ...
                     'LineWidth', 2, 'Color', 'red');
plot(dt(len-Size+2:end),widthchartdata1fhs(:,8), 'LineStyle', '-.', ...
                     'LineWidth', 2, 'Color', 'blue');              
legend('Sp500 Ret','VaR Interval interval','ES Interval Width');
legend('boxoff')
set(gcf, 'color', 'white');

subplot(2,2,3)
TitleStr = 'SP500 returns with 1% iid-Bootstrap tGJR VaR and ES interval width';
title(TitleStr) 
grid on;
hold on;

plot(dt(len-Size+2:end),iidwidthchartdata1(:,1), 'LineStyle', '--' , ...
                      'LineWidth', 1, 'Color', 'cyan');
plot(dt(len-Size+2:end),iidwidthchartdata1(:,5), 'LineStyle', '-.', ...
                     'LineWidth', 2, 'Color', 'red');
plot(dt(len-Size+2:end),iidwidthchartdata1(:,9), 'LineStyle', '-.', ...
                     'LineWidth', 2, 'Color', 'blue');     
legend('Sp500 Ret','VaR Interval interval','ES Interval Width');
legend('boxoff')
set(gcf, 'color', 'white');


%FHS plot
subplot(2,2,4)
TitleStr = 'SP500 returns with 1% iid-FHS tGJR VaR and ES interval width';
title(TitleStr) 
grid on;
hold on;
plot(dt(len-Size+2:end),widthchartdata1fhs(:,1), 'LineStyle', '--' , ...
                      'LineWidth', 1, 'Color', 'cyan');
plot(dt(len-Size+2:end),widthchartdata1fhs(:,5), 'LineStyle', '-.', ...
                     'LineWidth', 2, 'Color', 'red');
plot(dt(len-Size+2:end),widthchartdata1fhs(:,9), 'LineStyle', '-.', ...
                     'LineWidth', 2, 'Color', 'blue');              
legend('Sp500 Ret','VaR Interval interval','ES Interval Width');
legend('boxoff')
set(gcf, 'color', 'white');



% 5%
figure(2)
subplot(2,2,1)
TitleStr = 'SP500 returns with 5% iid-Bootstrap GARCH VaR and ES interval width';
title(TitleStr) 
grid on;
hold on;

plot(dt(len-Size+2:end),iidwidthchartdata5(:,1), 'LineStyle', '--' , ...
                      'LineWidth', 1, 'Color', 'cyan');
plot(dt(len-Size+2:end),iidwidthchartdata5(:,2), 'LineStyle', '-.', ...
                     'LineWidth', 2, 'Color', 'red');
plot(dt(len-Size+2:end),iidwidthchartdata5(:,6), 'LineStyle', '-.', ...
                     'LineWidth', 2, 'Color', 'blue');     
legend('Sp500 Ret','VaR Interval interval','ES Interval Width');
legend('boxoff')
set(gcf, 'color', 'white');


%FHS plot
subplot(2,2,2)
TitleStr = 'SP500 returns with 5% iid-FHS GARCH VaR and ES interval width';
title(TitleStr) 
grid on;
hold on;
plot(dt(len-Size+2:end),widthchartdata5fhs(:,1), 'LineStyle', '--' , ...
                      'LineWidth', 1, 'Color', 'cyan');
plot(dt(len-Size+2:end),widthchartdata5fhs(:,2), 'LineStyle', '-.', ...
                     'LineWidth', 2, 'Color', 'red');
plot(dt(len-Size+2:end),widthchartdata5fhs(:,6), 'LineStyle', '-.', ...
                     'LineWidth', 2, 'Color', 'blue');              
legend('Sp500 Ret','VaR Interval interval','ES Interval Width');
legend('boxoff')
set(gcf, 'color', 'white');

subplot(2,2,3)
TitleStr = 'SP500 returns with 5% iid-Bootstrap GJR VaR and ES interval width';
title(TitleStr) 
grid on;
hold on;

plot(dt(len-Size+2:end),iidwidthchartdata5(:,1), 'LineStyle', '--' , ...
                      'LineWidth', 1, 'Color', 'cyan');
plot(dt(len-Size+2:end),iidwidthchartdata5(:,3), 'LineStyle', '-.', ...
                     'LineWidth', 2, 'Color', 'red');
plot(dt(len-Size+2:end),iidwidthchartdata5(:,7), 'LineStyle', '-.', ...
                     'LineWidth', 2, 'Color', 'blue');     
legend('Sp500 Ret','VaR Interval interval','ES Interval Width');
legend('boxoff')
set(gcf, 'color', 'white');


%FHS plot
subplot(2,2,4)
TitleStr = 'SP500 returns with 5% iid-FHS GJR VaR and ES interval width';
title(TitleStr) 
grid on;
hold on;
plot(dt(len-Size+2:end),widthchartdata5fhs(:,1), 'LineStyle', '--' , ...
                      'LineWidth', 1, 'Color', 'cyan');
plot(dt(len-Size+2:end),widthchartdata5fhs(:,3), 'LineStyle', '-.', ...
                     'LineWidth', 2, 'Color', 'red');
plot(dt(len-Size+2:end),widthchartdata5fhs(:,7), 'LineStyle', '-.', ...
                     'LineWidth', 2, 'Color', 'blue');              
legend('Sp500 Ret','VaR Interval interval','ES Interval Width');
legend('boxoff')
set(gcf, 'color', 'white');



%student-t
figure(3)
subplot(2,2,1)
TitleStr = 'SP500 returns with 5% iid-Bootstrap tGARCH VaR and ES interval width';
title(TitleStr) 
grid on;
hold on;

plot(dt(len-Size+2:end),iidwidthchartdata5(:,1), 'LineStyle', '--' , ...
                      'LineWidth', 1, 'Color', 'cyan');
plot(dt(len-Size+2:end),iidwidthchartdata5(:,4), 'LineStyle', '-.', ...
                     'LineWidth', 2, 'Color', 'red');
plot(dt(len-Size+2:end),iidwidthchartdata5(:,8), 'LineStyle', '-.', ...
                     'LineWidth', 2, 'Color', 'blue');     
legend('Sp500 Ret','VaR Interval interval','ES Interval Width');
legend('boxoff')
set(gcf, 'color', 'white');


%FHS plot
subplot(2,2,2)
TitleStr = 'SP500 returns with 5% iid-FHS tGARCH VaR and ES interval width';
title(TitleStr) 
grid on;
hold on;
plot(dt(len-Size+2:end),widthchartdata5fhs(:,1), 'LineStyle', '--' , ...
                      'LineWidth', 1, 'Color', 'cyan');
plot(dt(len-Size+2:end),widthchartdata5fhs(:,4), 'LineStyle', '-.', ...
                     'LineWidth', 2, 'Color', 'red');
plot(dt(len-Size+2:end),widthchartdata5fhs(:,8), 'LineStyle', '-.', ...
                     'LineWidth', 2, 'Color', 'blue');              
legend('Sp500 Ret','VaR Interval interval','ES Interval Width');
legend('boxoff')
set(gcf, 'color', 'white');

subplot(2,2,3)
TitleStr = 'SP500 returns with 5% iid-Bootstrap tGJR VaR and ES interval width';
title(TitleStr) 
grid on;
hold on;

plot(dt(len-Size+2:end),iidwidthchartdata5(:,1), 'LineStyle', '--' , ...
                      'LineWidth', 1, 'Color', 'cyan');
plot(dt(len-Size+2:end),iidwidthchartdata5(:,5), 'LineStyle', '-.', ...
                     'LineWidth', 2, 'Color', 'red');
plot(dt(len-Size+2:end),iidwidthchartdata5(:,9), 'LineStyle', '-.', ...
                     'LineWidth', 2, 'Color', 'blue');     
legend('Sp500 Ret','VaR Interval interval','ES Interval Width');
legend('boxoff')
set(gcf, 'color', 'white');


%FHS plot
subplot(2,2,4)
TitleStr = 'SP500 returns with 5% iid-FHS tGJR- VaR and ES interval width';
title(TitleStr) 
grid on;
hold on;
plot(dt(len-Size+2:end),widthchartdata5fhs(:,1), 'LineStyle', '--' , ...
                      'LineWidth', 1, 'Color', 'cyan');
plot(dt(len-Size+2:end),widthchartdata5fhs(:,5), 'LineStyle', '-.', ...
                     'LineWidth', 2, 'Color', 'red');
plot(dt(len-Size+2:end),widthchartdata5fhs(:,9), 'LineStyle', '-.', ...
                     'LineWidth', 2, 'Color', 'blue');              
legend('Sp500 Ret','VaR Interval interval','ES Interval Width');
legend('boxoff')
set(gcf, 'color', 'white');



