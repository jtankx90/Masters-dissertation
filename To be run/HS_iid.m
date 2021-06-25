clear all

oridatasp500 = readtable('dissertation data-finalv2.xlsx','Sheet','SP500');
ret = [table2array(oridatasp500(:,3))*100];

len = [length(ret)];

%Note: require manual change of backtesting window size to 250 and 1250
%samplesize=len-1250
samplesize= len-250;
Size = len-samplesize;
index = [0.05,0.01];




B = 1000;
cutoff = round(index.*samplesize);
VaR5 = NaN(B,Size);
ES5 = NaN(B,Size);

VaR1 = NaN(B,Size);
ES1 = NaN(B,Size);

for i=1:Size
   display(i) 
   window = ret(i:(samplesize+i-1));

   parfor k = 1:B
       bootwindow = datasample(window,samplesize);
       sortbootwindow = sort(bootwindow)
       VaR5(k,i) = -sortbootwindow(cutoff(1));
       ES5(k,i) = -mean(sortbootwindow(1:cutoff(1)));

       VaR1(k,i) = -sortbootwindow(cutoff(2));
       ES1(k,i) = -mean(sortbootwindow(1:cutoff(2)));

   end
   
end

avgVaR = [mean(VaR5,1).',mean(VaR1,1).'];
avgES = [mean(ES5,1).',mean(ES1,1).'];


stdVaR = [std(VaR5,1).',std(VaR1,1).'];
stdES = [std(ES5,1).',std(ES1,1).'];

avgSTDvar = mean(stdVaR);
avgSTDes = mean(stdES);

%% Prediction Interval
PI = 0.05;
bound = [round(PI/2*B),round((1-PI/2)*B)];

%VaR
VaRPI5 = sort(VaR5,1);
VaRPI5 = VaRPI5(bound,:).';

VaRPI1 = sort(VaR1,1);
VaRPI1 = VaRPI1(bound,:).';

combinedVaRPI = [VaRPI5, VaRPI1];

%ES
ESPI5 = sort(ES5,1);
ESPI5 = ESPI5(bound,:).';

ESPI1 = sort(ES1,1);
ESPI1 = ESPI1(bound,:).';

combinedESPI = [ESPI5, ESPI1];


%Average PI
avgVarPI5 = mean(VaRPI5,1);
avgVarPI1 = mean(VaRPI1,1);
avgESPI5 = mean(ESPI5,1);
avgESPI1 = mean(ESPI1,1);


%Width of Interval
%VaR
varPIlen5 = minus(VaRPI5(:,2),VaRPI5(:,1));
varPIlen1 = minus(VaRPI1(:,2),VaRPI1(:,1));


%ES
esPIlen5 = minus(ESPI5(:,2),ESPI5(:,1));
esPIlen1 = minus(ESPI1(:,2),ESPI1(:,1));




%Standard Error
stdVaR = [std(VaR5,1).',std(VaR1,1).'];
stdES = [std(ES5,1).',std(ES1,1).'];

avgSTDvar = mean(stdVaR);
avgSTDes = mean(stdES);






%% Back-testing
Breach = NaN(Size,2);

for i=1:2
    Breach(:,i) = ret(samplesize+1:end) < -avgVaR(:,i);
end
vio5 = sum(Breach(:,1));
vio1 = sum(Breach(:,2));
vr5= vio5./(index(1)*Size);
vr1= vio1./(index(2)*Size);



% Unconditional Coverage Test
%5 percent
phat(:,1) = vio5/Size;
phat(:,2) = vio1/Size;
V1(:,1) = vio5;
V1(:,2) = vio1;
V0= Size - V1;

for i =1:2
    LogLConst(:,i) = V1(:,i)*log(index(i)) + V0(:,i)*log(1-index(i));
    if phat(:,i) == 0
         LogLUnconst(:,i) = V0(:,i)*log(1-phat(:,i));
    else
         LogLUnconst(:,i) = V1(:,i)*log(phat(:,i)) + V0(:,i)*log(1-phat(:,i));
    end
end
for i = 1:2
    TestStatisticUC(1,i)          = -2 * (LogLConst(:,i) - LogLUnconst(:,i));
    SignificanceUC(:,i)      = 1 - chi2cdf(TestStatisticUC(1,i), 1);
end

% Independence Test
ind = NaN(1,2);
testsigcc = NaN(1,2);
ind(1,:) = [ind_test(Breach(:,1)),ind_test(Breach(:,2))];
testsigcc(1,:) = 1-chi2cdf(ind(1,:),1);


%Conditional Coverage
LRCC = NaN(1,2);
LRpvalues= NaN(1,2);
for i =1:2;
   LRCC(1,i) =   TestStatisticUC(1,i) + ind(1,i);
   LRpvalues(1,i) = 1-chi2cdf(LRCC(1,i),2);
end
%% Backtesting ES
retwindow = ret(samplesize+1:end);

for i =1:2
    q = find(retwindow <= -avgVaR(:,i));
    nES(:,i) = mean(retwindow(q) ./ -avgES(q,i));
end


%% Plot
dt = datetime(table2array(oridatasp500(:,1)),'InputFormat','dd/MM/yyyy');



Figno = 1
% 5 percent
figure(Figno);
subplot(2,1,1)
chartdata = [retwindow,-avgVaR(:,1),-avgES(:,1)];
TitleStr = strcat({'S&P Returns with '},{'5% HS VaR and ES '});
title(TitleStr) 
grid on;
hold on;
plot(dt(len-Size+1:end),chartdata(:,1), 'LineStyle', '--' , ...
                  'LineWidth', 0.5, 'Color', 'black');
plot(dt(len-Size+1:end),chartdata(:,2), 'LineStyle', '-.', ...
                 'LineWidth', 1, 'Color', 'red');
plot(dt(len-Size+1:end),chartdata(:,3), 'LineStyle', '-', ...
                 'LineWidth', 1, 'Color', 'Blue');   

legend('S&P Returns','VaR','ES');
legend('boxoff')
set(gcf, 'color', 'white');
    




subplot(2,1,2)
chartdata = [retwindow,-avgVaR(:,2),-avgES(:,2)];
TitleStr = strcat({'S&P Returns with '},{'1% HS VaR and ES '});
title(TitleStr) 
grid on;
hold on;
plot(dt(len-Size+1:end),chartdata(:,1), 'LineStyle', '--' , ...
                  'LineWidth', 0.5, 'Color', 'black');
plot(dt(len-Size+1:end),chartdata(:,2), 'LineStyle', '-.', ...
                 'LineWidth', 1, 'Color', 'red');
plot(dt(len-Size+1:end),chartdata(:,3), 'LineStyle', '-', ...
                 'LineWidth', 1, 'Color', 'Blue');   

legend('S&P Returns','VaR','ES');
legend('boxoff')
set(gcf, 'color', 'white');


   
% Width difference plot
HSlenchartdata1 = [retwindow,(minus(esPIlen1(:,:),varPIlen1(:,:)))];
HSlenchartdata5 = [retwindow,(minus(esPIlen5(:,:),varPIlen5(:,:)))];

dt = datetime(table2array(oridatasp500(:,1)),'InputFormat','dd/MM/yyyy');

figure(1)
subplot(2,1,1)
TitleStr = 'SP500 returns with width difference between 1% iid-HS VaR and ES';
title(TitleStr) 
grid on;
hold on;
plot(dt(len-Size+1:end),HSlenchartdata1(:,1), 'LineStyle', '-' , ...
                      'LineWidth', 1, 'Color', 'cyan');
plot(dt(len-Size+1:end),HSlenchartdata1(:,2), 'LineStyle', '-.', ...
                     'LineWidth', 2, 'Color', 'red');

legend('Sp500 Ret','Width Diff');
legend('boxoff')
set(gcf, 'color', 'white');

subplot(2,1,2)
TitleStr = 'SP500 returns with width difference between 5% iid-HS VaR and ES';
title(TitleStr) 
grid on;
hold on;
plot(dt(len-Size+1:end),HSlenchartdata5(:,1), 'LineStyle', '-' , ...
                      'LineWidth', 1, 'Color', 'cyan');
plot(dt(len-Size+1:end),HSlenchartdata5(:,2), 'LineStyle', '-.', ...
                     'LineWidth', 2, 'Color', 'red');           
legend('Sp500 Ret','Width Diff');
legend('boxoff')
set(gcf, 'color', 'white');


% Width Plot

HSwidthchartdata1 = [retwindow,esPIlen1(:,:),varPIlen1(:,:)];
HSwidthchartdata5 = [retwindow,esPIlen5(:,:),varPIlen5(:,:)];


figure(2)
subplot(2,1,1)
TitleStr = 'SP500 returns with width of 1% iid-HS VaR and ES';
title(TitleStr) 
grid on;
hold on;
plot(dt(len-Size+1:end),HSwidthchartdata1(:,1), 'LineStyle', '-' , ...
                      'LineWidth', 1, 'Color', 'cyan');
plot(dt(len-Size+1:end),HSwidthchartdata1(:,2), 'LineStyle', '-.', ...
                     'LineWidth', 2, 'Color', 'blue');
plot(dt(len-Size+1:end),HSwidthchartdata1(:,3), 'LineStyle', '-.', ...
                     'LineWidth', 2, 'Color', 'red');

legend('Sp500 Ret','ES','VaR');
legend('boxoff')
set(gcf, 'color', 'white');

subplot(2,1,2)
TitleStr = 'SP500 returns with width of 5% iid-HS VaR and ES';
title(TitleStr) 
grid on;
hold on;
plot(dt(len-Size+1:end),HSwidthchartdata5(:,1), 'LineStyle', '-' , ...
                      'LineWidth', 1, 'Color', 'cyan');
plot(dt(len-Size+1:end),HSwidthchartdata5(:,2), 'LineStyle', '-.', ...
                     'LineWidth', 2, 'Color', 'blue');
plot(dt(len-Size+1:end),HSwidthchartdata5(:,3), 'LineStyle', '-.', ...
                     'LineWidth', 2, 'Color', 'red');          
legend('Sp500 Ret','ES','VaR');
legend('boxoff')
set(gcf, 'color', 'white');
