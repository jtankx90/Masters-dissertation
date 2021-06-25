clear all

oridatasp500 = readtable('dissertation data-finalv2.xlsx','Sheet','SP500');
ret = [table2array(oridatasp500(:,3))*100];
date = [table2array(oridatasp500(:,1))];
len = [length(ret)];
date2 = datetime(date,'InputFormat','dd/MM/yyyy');

indsize = [250,500,750,1250,2000,2500];

win = len-indsize;
%Manual change to vary estimation window size
% samplesize= len-win(1);
% samplesize= len-win(2);
% samplesize= len-win(3);
% samplesize= len-win(4);
% samplesize= len-win(5);
samplesize= len-win(6);

Size = len-samplesize;
index = [0.05,0.01];

cutoff = round(index.*samplesize);
VaR = NaN(Size,2);
ES = NaN(Size,2);

for i=1:Size
   window(:,i) = ret(i:(samplesize+i-1));
   sortwindow = sort(window(:,i));
   VaR(i,1) = -sortwindow(cutoff(1));
   ES(i,1) = -mean(sortwindow(1:cutoff(1)));
    
   VaR(i,2) = -sortwindow(cutoff(2));
   ES(i,2) = -mean(sortwindow(1:cutoff(2)));
    
    
end
combinedchart = [-VaR,-ES];

%% Back-testing
Breach = NaN(Size,2);

for i=1:2
    Breach(:,i) = ret(samplesize+1:end) < -VaR(:,i);
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
    TestStatisticUC(i,1)          = -2 * (LogLConst(:,i) - LogLUnconst(:,i));
    SignificanceUC(:,i)      = 1 - chi2cdf(TestStatisticUC(i,1), 1);
end


% Independence Test

ind = NaN(2,1);
testsigcc = NaN(1,2);
for i =1:2;
    ind(i,1) = ind_test(Breach(:,i));
    testsigcc(1,i) = 1-chi2cdf(ind(i,1),1);
end;



%CC test

LRCC = NaN(2,1);
LRpvalues= NaN(1,2);

for i =1:2;
   LRCC(i,1) =   TestStatisticUC(i,1) + ind(i,1);
   LRpvalues(1,i) = 1-chi2cdf(LRCC(i,1),2);
end

%% Backtesting ES
retwindow = ret(samplesize+1:end);


for i =1:2
    q = find(retwindow <= -VaR(:,i));
    nES(:,i) = mean(retwindow(q) ./ -ES(q,i));
end










% load data for plot
ten = load('chartdatawindow2500.mat');
ten = ten.combinedchart;
tenvar = ten(:,1:2);
tenes = ten(:,3:4);

eight = load('chartdatawindow2000.mat');
eight = eight.combinedchart;
eightvar = eight(:,1:2);
eightes = eight(:,3:4);


five = load('chartdatawindow1250.mat');
five =five.combinedchart;
fivevar = five(:,1:2);
fivees = five(:,3:4);



three = load('chartdatawindow750.mat');
three = three.combinedchart;
threevar = three(:,1:2);
threes = three(:,3:4);


two = load('chartdatawindow500.mat');
two = two.combinedchart;
twovar = two(:,1:2);
twoes = two(:,3:4);

on = load('chartdatawindow250.mat');
on = on.combinedchart;
onvar =on(:,1:2);
oness = on(:,3:4);

% HS plot
alpha = {'5%','1%'}
Figno=1
for i=1:3
   figure(Figno);
    axes('Units', 'normalized', 'Position', [0 0 1 1])
    if i==1
      for j = 1:2
            retwindow10 = ret(len-length(ten)+1:end);
            date10 = date2(len-length(ten)+1:end);
            subplot(3,2,[j])
            chartdata10 = [retwindow10,tenvar(:,j),tenes(:,j)];
            TitleStr = strcat({'S&P Returns with '},alpha(j),' HS-VaR ', ' and HS-ES');
            title(TitleStr) 
            grid on;
            hold on;
            plot(date10,chartdata10(:,1), 'LineStyle', '--' , ...
                              'LineWidth', 0.5, 'Color', 'black');
            plot(date10,chartdata10(:,2), 'LineStyle', '-.', ...
                             'LineWidth', 1, 'Color', 'red');
            plot(date10,chartdata10(:,3), 'LineStyle', '-.', ...
                             'LineWidth', 1, 'Color', 'blue');
            legend(['S&P Returns',strcat(alpha(j), ' HS-VaR'),strcat(alpha(j),' HS-ES')],'Location','best');
            legend('boxoff')  
            if j ==2
                ylabel('')
            else
                ylabel('Samplesize = 2500')
            end
     

          
            set(gcf, 'color', 'white');
      end
  
    end
    
    if i ==2
        for j = 1:2
            retwindow8 = ret((len-length(eight)+1):end);
            date8 = date2(len-length(eight)+1:end);
            subplot(3,2,[i+j])
            chartdata8 = [retwindow8,eightvar(:,j),eightes(:,j)];
            grid on;
            hold on;
            plot(date8,chartdata8(:,1), 'LineStyle', '--' , ...
                              'LineWidth', 0.5, 'Color', 'black');
            plot(date8,chartdata8(:,2), 'LineStyle', '-.', ...
                             'LineWidth', 1, 'Color', 'red');
            plot(date8,chartdata8(:,3), 'LineStyle', '-.', ...
                             'LineWidth', 1, 'Color', 'blue');
              if j ==2
                ylabel('')
            else
                ylabel('Samplesize = 2000')
            end

            set(gcf, 'color', 'white');
        end
    end
    
      if i ==3
        for j = 1:2
            retwindow5 = ret(len-length(five)+1:end);
            date5 = date2(len-length(five)+1:end);
            subplot(3,2,[i+j+1])
            chartdata5 = [retwindow5,fivevar(:,j),fivees(:,j)];
            grid on;
            hold on;
            plot(date5,chartdata5(:,1), 'LineStyle', '--' , ...
                              'LineWidth', 0.5, 'Color', 'black');
            plot(date5,chartdata5(:,2), 'LineStyle', '-.', ...
                             'LineWidth', 1, 'Color', 'red');
            plot(date5,chartdata5(:,3), 'LineStyle', '-.', ...
                             'LineWidth', 1, 'Color', 'blue');
              if j ==2
                ylabel('')
            else
                ylabel('Samplesize = 1250')
            end

           
            set(gcf, 'color', 'white');
        end
      end
     


end


for i=1:3
   figure(2);
    axes('Units', 'normalized', 'Position', [0 0 1 1])
  if i ==1
        for j = 1:2
            
            retwindow3 = ret(len-length(three)+1:end);
            date3 = date2(len-length(three)+1:end);
            subplot(3,2,[j])
            chartdata3 = [retwindow3,threevar(:,j),threes(:,j)];
              TitleStr = strcat({'S&P Returns with '},alpha(j),' HS-VaR ', ' and HS-ES');
            title(TitleStr) 
            grid on;
            hold on;
            plot(date3,chartdata3(:,1), 'LineStyle', '--' , ...
                              'LineWidth', 0.5, 'Color', 'black');
            plot(date3,chartdata3(:,2), 'LineStyle', '-.', ...
                             'LineWidth', 1, 'Color', 'red');
            plot(date3,chartdata3(:,3), 'LineStyle', '-.', ...
                             'LineWidth', 1, 'Color', 'blue');
            legend(['S&P Returns',strcat(alpha(j), ' HS-VaR'),strcat(alpha(j),' HS-ES')],'Location','best');
            legend('boxoff')
              if j ==2
                ylabel('')
            else
                ylabel('Samplesize = 750')
            end

      
            set(gcf, 'color', 'white');
        end
    end
      if i ==2
        for j = 1:2
            retwindow2 = ret(len-length(two)+1:end);
            datewin2 = date2(len-length(two)+1:end);
            subplot(3,2,[i+j])
            chartdata2 = [retwindow2,twovar(:,j),twoes(:,j)];
            grid on;
            hold on;
            plot(datewin2,chartdata2(:,1), 'LineStyle', '--' , ...
                              'LineWidth', 0.5, 'Color', 'black');
            plot(datewin2,chartdata2(:,2), 'LineStyle', '-.', ...
                             'LineWidth', 1, 'Color', 'red');
            plot(datewin2,chartdata2(:,3), 'LineStyle', '-.', ...
                             'LineWidth', 1, 'Color', 'blue');
            if j ==2
                ylabel('')
            else
                ylabel('Samplesize = 500')
            end

 
            set(gcf, 'color', 'white');
        end
     end
    
      if i ==3
        for j = 1:2
            retwindow = ret(len-length(on)+1:end);
            datewin = date2(len-length(on)+1:end);
            subplot(3,2,[i+j+1])
            chartdata = [retwindow,onvar(:,j),oness(:,j)];
            grid on;
            hold on;
            plot(datewi n,chartdata(:,1), 'LineStyle', '--' , ...
                              'LineWidth', 0.5, 'Color', 'black');
            plot(datewin,chartdata(:,2), 'LineStyle', '-.', ...
                             'LineWidth', 1, 'Color', 'red');
            plot(datewin,chartdata(:,3), 'LineStyle', '-.', ...
                             'LineWidth', 1, 'Color', 'blue');
                 if j ==2
                ylabel('')
            else
                ylabel('Samplesize = 250')
            end

    
            set(gcf, 'color', 'white'); 
        end
    end
end


