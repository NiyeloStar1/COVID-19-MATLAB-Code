%% Run the ODE model
%Keep running tidy
% clear all
% close all
clc

%% Set up the plot
cmap = colormap(parula(10));
set(gca,'fontsize',16)
set(0,'defaultaxeslinestyleorder',{'-*',':','-','.'}) %or whatever you want
set(0,'defaultaxesfontsize',16)
set(0,'defaultlinelinewidth',2)

%% Load all data needed at this point

Polymod = load("MixingData.mat");
varyingpara = load("parameters.mat");
prob = load("Probab.mat");

% Load COVID-19 test data
Tdata = readtable('test_data.xlsx','Sheet','Jan_SepTest');

%% Load the population data 

Pop= load("EnglandDecilePop.mat");

N =  Pop.TotalPop;
D = (sum(N,2))';        %Total by socioeconomic group
N_age = (sum(N,1));     %Total by Age group
total = sum(D);         %Total England Population

N1= N(1,:); N2= N(2,:);N3= N(3,:);N4= N(4,:);N5= N(5,:);
N6= N(6,:); N7= N(7,:);N8= N(8,:);N9= N(9,:);N10= N(10,:);

%% Extract the age mixing matrix

Mix_H = Polymod.UK_from_toH;
Mix_S = Polymod.UK_from_toS;
Mix_W = Polymod.UK_from_toW;
Mix_O = Polymod.UK_from_toO;

%% Extract the infection parameter by Age

a = varyingpara.Susceptibility; % Susceptibility
h = varyingpara.Detection; % Infectivity

%% Setup the social mixing matrix
x = [10 9 8 7 6 5 4 3 2 1];


epi = 0.3;
soc = SocioContMix(x,epi,total,D);

soc = soc/sum(sum(soc));
%% Define the number of age and social groups

% soc = 10;   %Sociodemographic group
n = 21;     %age group

%% Set up the Parameter

para = struct('tau',1, ...
    'epsilon',2/100, ...
    'p',0.8526, ...
    'rho',0.4,'theta',0.6,...
    'delta',prob.Sympt_2_hosp, ...
    'gamma',1/14,'gamma_2',1/14, ...
    'd',prob.Hosp_2_Death, ...
    'PHI',3,...
    'a',a,'h',h,'N',Polymod.UK_PP,...
    'beta',0.099,...%0.099
    'iota',0.27,...
    'N_total',total,'N1',N1,'N2',N2,'N3',N3,'N4',N4, ...
    'N5',N5,'N6',N6,'N7',N7,'N8',N8,'N9',N9,'N10',N10,'N_age',N_age,'D',N, ...
    'N_social',D);


%%%%%%%% CALCULATE R0
% ageMix = Mix_O+Mix_W+Mix_S+Mix_H;


% for i =1:10
% R0(i,:) = sum(para.beta*((sum(soc(i,:))*(para.a*ageMix.*para.h)))*1/21);
% end
% R0_total = sum(R0,2);
% 


%% Set up the run-time
maxtime = 200;

%% Define the initial conditions
EICs = load("ExposedIcs.mat"); %Load the exposed matrix
%Extract the exposed class
E = EICs.Exposed; %Exposed
I1 =[0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0];
I2 = E(2,:);I3 =E(3,:);I4 =E(4,:);I5 =E(5,:);
I6 =E(6,:);I7 =E(7,:);I8 =E(8,:);I9 =E(9,:);I10 =E(10,:);
% I1 = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
% J1 = [0 0 0 0 0 0 0 0 0 0 0 0 5 0 0 0 0 0 0 0 0];

% Susceptible
S1 = N1 -I1;
% S2 = N2 -E2; S3 = N3 -E3; S4 = N4 -E4; S5 = N5 -E5; 
% S6 = N6 -E6;S7 = N7 -E7; S8 = N8 -E8; S9 = N9 -E9; S10 = N10 -E10; 



ICs = struct('S1',S1,'E1',zeros(1,n),'A1',zeros(1,n),'J1',zeros(1,n), ...
    'Q1',zeros(1,n),'I1',I1,'H1',zeros(1,n),'R1',zeros(1,n), ...
    'D1',zeros(1,n),'Ir1',zeros(1,n),'S2',N2,'E2',zeros(1,n),'A2',zeros(1,n),'J2',zeros(1,n), ...
    'Q2',zeros(1,n),'I2',I2,'H2',zeros(1,n),'R2',zeros(1,n), ...
    'D2',zeros(1,n),'Ir2',zeros(1,n),'S3',N3,'E3',zeros(1,n),'A3',zeros(1,n),'J3',zeros(1,n), ...
    'Q3',zeros(1,n),'I3',I3,'H3',zeros(1,n),'R3',zeros(1,n), ...
    'D3',zeros(1,n),'Ir3',zeros(1,n),'S4',N4,'E4',zeros(1,n),'A4',zeros(1,n),'J4',zeros(1,n), ...
    'Q4',zeros(1,n),'I4',I4,'H4',zeros(1,n),'R4',zeros(1,n), ...
    'D4',zeros(1,n),'Ir4',zeros(1,n),'S5',N5,'E5',zeros(1,n),'A5',zeros(1,n),'J5',zeros(1,n), ...
    'Q5',zeros(1,n),'I5',I5,'H5',zeros(1,n),'R5',zeros(1,n), ...
    'D5',zeros(1,n),'Ir5',zeros(1,n),'S6',N6,'E6',zeros(1,n),'A6',zeros(1,n),'J6',zeros(1,n), ...
    'Q6',zeros(1,n),'I6',I6,'H6',zeros(1,n),'R6',zeros(1,n), ...
    'D6',zeros(1,n),'Ir6',zeros(1,n),'S7',N7,'E7',zeros(1,n),'A7',zeros(1,n),'J7',zeros(1,n), ...
    'Q7',zeros(1,n),'I7',I7,'H7',zeros(1,n),'R7',zeros(1,n), ...
    'D7',zeros(1,n),'Ir7',zeros(1,n),'S8',N8,'E8',zeros(1,n),'A8',zeros(1,n),'J8',zeros(1,n), ...
    'Q8',zeros(1,n),'I8',I8,'H8',zeros(1,n),'R8',zeros(1,n), ...
    'D8',zeros(1,n),'Ir8',zeros(1,n),'S9',N9,'E9',zeros(1,n),'A9',zeros(1,n),'J9',zeros(1,n), ...
    'Q9',zeros(1,n),'I9',I9,'H9',zeros(1,n),'R9',zeros(1,n), ...
    'D9',zeros(1,n),'Ir9',zeros(1,n),'S10',N10,'E10',zeros(1,n),'A10',zeros(1,n),'J10',zeros(1,n), ...
    'Q10',zeros(1,n),'I10',I10,'H10',zeros(1,n),'R10',zeros(1,n), ...
    'D10',zeros(1,n),'Ir10',zeros(1,n));



%% Run the ODE
tic
[Classes,R0,R02] = ODE_social1T(para,ICs,maxtime,Mix_H,Mix_W,Mix_S,Mix_O,n,soc,x);
toc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TOTAL PREVALENCE BY SOCIODEMOGRAPHIC GROUP %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% GROUP ONE

prev1 = sum((Classes.A1 + Classes.J1 + Classes.Q1 +Classes.I1 +Classes.H1),2);

%% GROUP TWO

prev2 = sum((Classes.A2 + Classes.J2 + Classes.Q2 +Classes.I2 +Classes.H2),2);

%% GROUP THREE

prev3 = sum((Classes.A3 + Classes.J3 + Classes.Q3 +Classes.I3 +Classes.H3),2);

%% GROUP FOUR

prev4 = sum((Classes.A4 + Classes.J4 + Classes.Q4 +Classes.I4 +Classes.H4),2);

%% GROUP FIVE

prev5 = sum((Classes.A5 + Classes.J5 + Classes.Q5 +Classes.I5 +Classes.H5),2);

%% GROUP SIX

prev6 = sum((Classes.A6 + Classes.J6 + Classes.Q6 +Classes.I6 +Classes.H6),2);

%% GROUP SEVEN

prev7 = sum((Classes.A7 + Classes.J7 + Classes.Q7 +Classes.I7 +Classes.H7),2);

%% GROUP EIGHT

prev8 = sum((Classes.A8 + Classes.J8 + Classes.Q8 +Classes.I8 +Classes.H8),2);

%% GROUP NINE

prev9 = sum((Classes.A9 + Classes.J9 + Classes.Q9 +Classes.I9 +Classes.H9),2);

%% GROUP TEN

prev10 = sum((Classes.A10 + Classes.J10 + Classes.Q10 +Classes.I10 +Classes.H10),2);

PREV = [prev1,prev2,prev3,prev4,prev5,prev6,prev7,prev8,prev9,prev10];

% % Extract the diagonal elements
% diagonal_points = diag(R0);
% % Plot the diagonal points
% 
% figure;
% hold on;
% for i = 1:numel(diagonal_points)
%     % Determine color based on the value of y
%     color_value = i / numel(diagonal_points); % Normalize the index
% 
%     % Determine circle size based on the value of y
%     circle_size = diagonal_points(i) * 50; % Adjust multiplier as needed
% 
%     % Plot the point
%     scatter(i, diagonal_points(i), circle_size, [color_value 0 1-color_value], 'filled');
% end
% 
% % Draw a horizontal dashed line at y=1
% plot([0 numel(diagonal_points)+1], [1 1], '--k');
% 
% hold off;
% xlabel('Deprivation Decile');
% ylabel('R0 Value');
% title('Reproduction number for each decile');
% xlim([1 numel(diagonal_points)]);
% xticks(1:numel(diagonal_points));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TOTAL INCIDENCE BY SOCIODEMOGRAPHIC GROUP %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IR = [sum(Classes.I1,2),sum(Classes.I2,2),sum(Classes.I3,2),sum(Classes.I4,2),...
    sum(Classes.I5,2),sum(Classes.I6,2),sum(Classes.I7,2),sum(Classes.I8,2),...
    sum(Classes.I9,2),sum(Classes.I10,2)];
% writematrix(IR,'infection.xlsx', 'Sheet', 'MixingTypes1')



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%% HOSPITALISED BY SOCIODEMOGRAPHIC GROUP %%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hosp = [sum(Classes.H1,2),sum(Classes.H2,2),sum(Classes.H3,2),sum(Classes.H4,2),...
    sum(Classes.H5,2),sum(Classes.H6,2),sum(Classes.H7,2),sum(Classes.H8,2),...
    sum(Classes.H9,2),sum(Classes.H10,2)];

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%% Death BY SOCIODEMOGRAPHIC GROUP %%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Death = [sum(Classes.D1,2),sum(Classes.D2,2),sum(Classes.D3,2),sum(Classes.D4,2),...
    sum(Classes.D5,2),sum(Classes.D6,2),sum(Classes.D7,2),sum(Classes.D8,2),...
    sum(Classes.D9,2),sum(Classes.D10,2)];

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%DAILY  NUMBER OF DEATH REPORTED PER DAY BY SOCIODEMOGRAPHIC GROUP %%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
%% Group One

De1 = [zeros(1,n);diff(Classes.D1)];
Total_De1 = sum(De1,2);

%% Group Two

De2 = [zeros(1,n);diff(Classes.D2)];
Total_De2 = sum(De2,2);

%% Group Three

De3 = [zeros(1,n);diff(Classes.D3)];
Total_De3 = sum(De3,2);

%% Group Four

De4 = [zeros(1,n);diff(Classes.D4)];
Total_De4 = sum(De4,2);

%% Group Five

De5 = [zeros(1,n);diff(Classes.D5)];
Total_De5 = sum(De5,2);

%% Group Six

De6 = [zeros(1,n);diff(Classes.D6)];
Total_De6 = sum(De6,2);

%% Group Seven

De7 = [zeros(1,n);diff(Classes.D7)];
Total_De7 = sum(De7,2);

%% Group Eight

De8 = [zeros(1,n);diff(Classes.D8)];
Total_De8 = sum(De8,2);

%% Group Nine

De9 = [zeros(1,n);diff(Classes.D9)];
Total_De9 = sum(De9,2);

%% Group Ten

De10 = [zeros(1,n);diff(Classes.D10)];
Total_De10 = sum(De10,2);

TOTAL_DE = [Total_De1,Total_De2,Total_De3,Total_De4,Total_De5,...
    Total_De6,Total_De7,Total_De8,Total_De9,Total_De10];



% writetable(table(TOTAL_DE), 'epi1.xlsx','sheet','Death')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%% COMBINE PLOT%%%%%%%%%%%%
xValues = 0:9;
figure;
subplot(2,2,1)
plot(Classes.t,PREV/10000)
% xticks(xValues);
% set(gca,'xtick',[]);
legend("1 - Most deprived",'2','3','4','5','6','7','8','9','10 - Least deprived')
% title('Total infection in each deprivation decile')
ylabel({'Total Infection'; '(per 10,000 people)'})
xlabel('Time(days)')
hold on
% xl = xline(84/30,'-r','DisplayName','Lowdown Period','LineWidth',3);
% xl = xline(164/30,'-g','DisplayName','Lowdown Easing','LineWidth',3);
hold off

subplot(2,2,2)
plot(Classes.t,IR/10000)
% xticks(xValues);
% xticklabels(customlabels);
% legend("1 - Most deprived",'2','3','4','5','6','7','8','9','10 - Least deprived')
ylabel({'Confirmed Cases'; '(per 10,000 people)'})
xlabel('Time(days)')
% title('Confirmed Cases in each deprivation decile')
hold on 
% xl = xline(84/30,'-r','DisplayName','Lowdown Period','LineWidth',3);
% xl = xline(164/30,'-g','DisplayName','Lowdown Easing','LineWidth',3);
hold off

subplot(2,2,3)
plot(Classes.t,Hosp/10000)
% xticks(xValues);
% set(gca,'xtick',[]);
% legend("1 - Most deprived",'2','3','4','5','6','7','8','9','10 - Least deprived')
% title('Hospitalised Individuals in each deprivation decile')
ylabel({'Hospitalised'; '(per 10,000 people)'})
xlabel('Time(days)')
hold on 
% xl = xline(84/30,'-r','DisplayName','Lowdown Period','LineWidth',3);
% xl = xline(164/30,'-g','DisplayName','Lowdown Easing','LineWidth',3);
hold off

subplot(2,2,4)
plot(Classes.t,TOTAL_DE/10000)
% xticks(xValues);
% legend("1 - Most deprived",'2','3','4','5','6','7','8','9','10 - Least deprived')
% title('Case fatality in each deprivation decile')
ylabel({'Daily Death'; '(per 10,000 people)'})
xlabel('Time(days)')
% xticklabels(customlabels);
hold on 
% xl = xline(84/30,'-r','DisplayName','Lowdown Period','LineWidth',3);
% xl = xline(164/30,'-g','DisplayName','Lowdown Easing','LineWidth',3);
hold off




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TOTAL PREVALENCE BY AGE GROUP %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% GROUP ONE

prevAge1 = ((Classes.A1 + Classes.J1 + Classes.Q1 +Classes.I1 +Classes.H1));

%% GROUP TWO

prevAge2 = ((Classes.A2 + Classes.J2 + Classes.Q2 +Classes.I2 +Classes.H2));

%% GROUP THREE

prevAge3 = ((Classes.A3 + Classes.J3 + Classes.Q3 +Classes.I3 +Classes.H3));

%% GROUP FOUR

prevAge4 = ((Classes.A4 + Classes.J4 + Classes.Q4 +Classes.I4 +Classes.H4));

%% GROUP FIVE

prevAge5 = ((Classes.A5 + Classes.J5 + Classes.Q5 +Classes.I5 +Classes.H5));

%% GROUP SIX

prevAge6 = ((Classes.A6 + Classes.J6 + Classes.Q6 +Classes.I6 +Classes.H6));

%% GROUP SEVEN

prevAge7 = ((Classes.A7 + Classes.J7 + Classes.Q7 +Classes.I7 +Classes.H7));

%% GROUP EIGHT

prevAge8 = ((Classes.A8 + Classes.J8 + Classes.Q8 +Classes.I8 +Classes.H8));

%% GROUP NINE

prevAge9 = ((Classes.A9 + Classes.J9 + Classes.Q9 +Classes.I9 +Classes.H9));

%% GROUP TEN

prevAge10 = ((Classes.A10 + Classes.J10 + Classes.Q10 +Classes.I10 +Classes.H10));

PREVAGE = prevAge1+prevAge2+prevAge3+prevAge4+prevAge5+prevAge6+prevAge7+prevAge8+prevAge9+prevAge10;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TOTAL INCIDENCE BY SOCIODEMOGRAPHIC GROUP %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IRAge = Classes.I1+Classes.I2+Classes.I3+Classes.I4+Classes.I5+Classes.I6+...
Classes.I7+Classes.I8+Classes.I9+Classes.I10;



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%% HOSPITALISED BY SOCIODEMOGRAPHIC GROUP %%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HospAge = Classes.H1+Classes.H2+Classes.H3+Classes.H4+Classes.H5+Classes.H6+...
Classes.H7+Classes.H8+Classes.H9+Classes.H10;


% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%% Death BY SOCIODEMOGRAPHIC GROUP %%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Death = Classes.D1+Classes.D2+Classes.D3+Classes.D4+Classes.D5+Classes.D6+...
Classes.D7+Classes.D8+Classes.D9+Classes.D10;

DeAGE = [zeros(1,n);diff(Death)];

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%% COMBINE PLOT%%%%%%%%%%%%

% xValues = 0:9;
% figure;
% subplot(2,2,1)
% plot(Classes.t/30,PREVAGE/10000)
% % xticks(xValues);
% % set(gca,'xtick',[]);
% legend('0-4', '5-9' ,'10-14', '15-19', '20-24', '25-29','30-34', ...
%         '35-39','40-44','45-49','50-54','55-59','60-64','65-69','70-74','75-79', ...
%         '80-84','85-89','90-94','95-99','99+', 'NumColumns', 2)
% % title('Total infection in each deprivation decile')
% ylabel({'Total Infection'; '(per 10,000 people)'})
% xlabel('Time(months)')
% hold on
% xl = xline(84/30,'-r','DisplayName','Lowdown Period','LineWidth',3);
% xl = xline(164/30,'-g','DisplayName','Lowdown Easing','LineWidth',3);
% hold off
% 
% subplot(2,2,2)
% plot(Classes.t/30,IRAge/10000)
% % xticks(xValues);
% % xticklabels(customlabels);
% % legend("1 - Most deprived",'2','3','4','5','6','7','8','9','10 - Least deprived')
% ylabel({'Confirmed Cases'; '(per 10,000 people)'})
% xlabel('Time(months)')
% % title('Confirmed Cases in each deprivation decile')
% hold on 
% xl = xline(84/30,'-r','DisplayName','Lowdown Period','LineWidth',3);
% xl = xline(164/30,'-g','DisplayName','Lowdown Easing','LineWidth',3);
% % legend('0-4', '5-9' ,'10-14', '15-19', '20-24', '25-29','30-34', ...
% %         '35-39','40-44','45-49','50-54','55-59','60-64','65-69','70-74','75-79', ...
% %         '80-84','85-89','90-94','95-99','99+', 'NumColumns', 2)
% hold off
% 
% subplot(2,2,3)
% plot(Classes.t/30,HospAge/10000)
% % xticks(xValues);
% % set(gca,'xtick',[]);
% % legend("1 - Most deprived",'2','3','4','5','6','7','8','9','10 - Least deprived')
% % title('Hospitalised Individuals in each deprivation decile')
% ylabel({'Hospitalised'; '(per 10,000 people)'})
% xlabel('Time(months)')
% hold on 
% xl = xline(84/30,'-r','DisplayName','Lowdown Period','LineWidth',3);
% xl = xline(164/30,'-g','DisplayName','Lowdown Easing','LineWidth',3);
% hold off
% 
% subplot(2,2,4)
% plot(Classes.t/30,DeAGE/10000)
% % xticks(xValues);
% % legend("1 - Most deprived",'2','3','4','5','6','7','8','9','10 - Least deprived')
% % title('Case fatality in each deprivation decile')
% ylabel({'Daily Death'; '(per 10,000 people)'})
% xlabel('Time(months)')
% % xticklabels(customlabels);
% hold on 
% xl = xline(84/30,'-r','DisplayName','Lowdown Period','LineWidth',3);
% xl = xline(164/30,'-g','DisplayName','Lowdown Easing','LineWidth',3);
% hold off



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the dot sizes based on the sorted values



% AA = max(PREV);
% AA2 = (sum(AA)/total)*100;
% PropAge = (AA./N_age)*100;
% 
% BB = max(DeAGE);
% ProDeath = (BB./N_age)*100
% 
% CC = max(HospAge);
% proHosp = (CC./N_age)*100
% 



% ix







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%

%Export Prevalence
% writetable(table(PREV), 'epi1.xlsx','sheet','PREV')

%% PLOT
% customlabels = {'Jan 3', 'Feb 3' ,'Mar 3', 'Apr 3', 'May 3', 'Jun 3','Jul 3','Aug 3','Sep 3'};
% xValues = 0:9;

% figure;
% plot(Classes.t/30,PREV)
% % xticks(xValues);
% % xticklabels(customlabels);
% xlabel('Time (months)')
% ylabel({'Total Number of Infectious'; '(per 100,000 people)'})
% legend("1",'2','3','4','5','6','7','8','9','10')
% hold on
% xl = xline(84/30,'-r','DisplayName','Lowdown Period','LineWidth',3);
% xl = xline(164/30,'-g','DisplayName','Lowdown Easing','LineWidth',3);
% hold off
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%
%Write to Table

% writetable(table(IR), 'epi1.xlsx','sheet','incidence')

% figure;
% plot(Classes.t/30,IR)
% % xticks(xValues);
% % xticklabels(customlabels);
% xlabel('Time (months)')
% ylabel({'Number of Reported Cases'; '(per 100,000 people)'})
% legend("1",'2','3','4','5','6','7','8','9','10')
% hold on 
% xl = xline(83/30,'-r','DisplayName','Lowdown Period','LineWidth',3);
% xl = xline(164/30,'-g','DisplayName','Lowdown Easing','LineWidth',3);
% hold off
% % 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure;
% subplot(2,1,1)
% plot(Classes.t/30,PREV/10000)
% % xticklabels(customlabels);
% % set(gca,'xtick',[]);
% legend("1 - Most deprived",'2','3','4','5','6','7','8','9','10 - Least deprived')
% % title('Total infection in each deprivation decile')
% ylabel({'Total Infection'; '(per 10,000 people)'})
% hold on
% xl = xline(84/30,'-r','DisplayName','Lowdown Period','LineWidth',3);
% xl = xline(164/30,'-g','DisplayName','Lowdown Easing','LineWidth',3);
% hold off
% 
% subplot(2,1,2)
% plot(Classes.t/30,IR/10000)
% % xticks(xValues);
% % xticklabels(customlabels);
% % legend("1 - Most deprived",'2','3','4','5','6','7','8','9','10 - Least deprived')
% ylabel({'Confirmed Cases'; '(per 10,000 people)'})
% xlabel('Time(months)')
% % title('Confirmed Cases in each deprivation decile')
% hold on 
% xl = xline(84/30,'-r','DisplayName','Lowdown Period','LineWidth',3);
% xl = xline(164/30,'-g','DisplayName','Lowdown Easing','LineWidth',3);
% hold off

% subplot(2,2,4)
% plot(Classes.t,IR/100000)
% xlabel('Time (days)')
% ylabel({'Cummulative Number of Reported Cases'; '(per 100,000 people)'})
% legend("1",'2','3','4','5','6','7','8','9','10')
% hold off
% 
% 
% 

%Write to table

% writetable(table(Hosp), 'epi1.xlsx','sheet','Hospital')

% figure(5)
% plot(Classes.t/30,Hosp)
% % xticks(xValues);
% % xticklabels(customlabels);
% xlabel('time (days)')
% ylabel({'Hospital Cases'; '(per 100,000 people)'})
% legend("1",'2','3','4','5','6','7','8','9','10')
% hold on 
% xl = xline(84/30,'-r','DisplayName','Lowdown Period','LineWidth',3);
% xl = xline(164/30,'-g','DisplayName','Lowdown Easing','LineWidth',3);
% hold off
% % 
% 



% % 
% % figure(6)
% % plot(Classes.t,Death/100000)
% % xlabel('time (days)')
% % ylabel({'Fatally Infectious'; '(per 100,000 people)'})
% % legend("1",'2','3','4','5','6','7','8','9','10')
% % hold off
% % 
% % 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FIND MAXIMUM VALUES

% % Find the maximum for each element across the 200 values
% IR_value = [sum(Classes.Ir1,2),sum(Classes.Ir2,2),sum(Classes.Ir3,2),sum(Classes.Ir4,2),...
%     sum(Classes.Ir5,2),sum(Classes.Ir6,2),sum(Classes.Ir7,2),sum(Classes.Ir8,2),...
%     sum(Classes.Ir9,2),sum(Classes.Ir10,2)];
%  max_values = (max(IR_value,[],1)./D)*100;
%  per_total= (sum(max(IR_value))/total)*100
% 
% 
% % Display the results
% for i = 1:numel(max_values)
%     fprintf('Maximum value for prev%d: %f\n', i, max_values(i));
% end
% 
% 
% 
% I = Classes.I1+Classes.I2+Classes.I3+Classes.I4+Classes.I5+Classes.I6+ ...
%     Classes.I7+Classes.I8+Classes.I9+Classes.I10;
% 
% H = Classes.H1+Classes.H2+Classes.H3+Classes.H4+Classes.H5+Classes.H6+ ...
%     Classes.H7+Classes.H8+Classes.H9+Classes.H10;
% 
% figure;
% plot(Classes.t/30,I/10000)
% ylabel({'Confirm Cases'; '(by age per 10,000 people)'})
% xlabel('Time(months)')
% legend('0-4', '5-9' ,'10-14', '15-19', '20-24', '25-29','30-34', ...
%         '35-39','40-44','45-49','50-54','55-59','60-64','65-69','70-74','75-79', ...
%         '80-84','85-89','90-94','95-99','99+')
% hold on 
% xl = xline(84/30,'-k','DisplayName','Lowdown Period','LineWidth',3);
% xl = xline(164/30,'-g','DisplayName','Lowdown Easing','LineWidth',3);
% hold off
% 




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %%%%%%DAILY  NUMBER OF DEATH REPORTED PER DAY BY SOCIODEMOGRAPHIC GROUP %%%%%%%%%%%%%%%%%%%%%%%%
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % 
% %% Group One
% 
% De1 = [zeros(1,n);diff(Classes.D1)];
% Total_De1 = sum(De1,2);
% 
% %% Group Two
% 
% De2 = [zeros(1,n);diff(Classes.D2)];
% Total_De2 = sum(De2,2);
% 
% %% Group Three
% 
% De3 = [zeros(1,n);diff(Classes.D3)];
% Total_De3 = sum(De3,2);
% 
% %% Group Four
% 
% De4 = [zeros(1,n);diff(Classes.D4)];
% Total_De4 = sum(De4,2);
% 
% %% Group Five
% 
% De5 = [zeros(1,n);diff(Classes.D5)];
% Total_De5 = sum(De5,2);
% 
% %% Group Six
% 
% De6 = [zeros(1,n);diff(Classes.D6)];
% Total_De6 = sum(De6,2);
% 
% %% Group Seven
% 
% De7 = [zeros(1,n);diff(Classes.D7)];
% Total_De7 = sum(De7,2);
% 
% %% Group Eight
% 
% De8 = [zeros(1,n);diff(Classes.D8)];
% Total_De8 = sum(De8,2);
% 
% %% Group Nine
% 
% De9 = [zeros(1,n);diff(Classes.D9)];
% Total_De9 = sum(De9,2);
% 
% %% Group Ten
% 
% De10 = [zeros(1,n);diff(Classes.D10)];
% Total_De10 = sum(De10,2);
% 
% TOTAL_DE = [Total_De1,Total_De2,Total_De3,Total_De4,Total_De5,...
%     Total_De6,Total_De7,Total_De8,Total_De9,Total_De10];
% 





% figure;
% plot(Classes.t,H)
% legend('0-4', '5-9' ,'10-14', '15-19', '20-24', '25-29','30-34', ...
%         '35-39','40-44','45-49','50-54','55-59','60-64','65-69','70-74','75-79', ...
%         '80-84','85-89','90-94','95-99','99+')
% 
% 
% 



















% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%%%%% HEAT MAPS %%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Hospitalised
% 
Hospital = [sum(Classes.H1,1);sum(Classes.H2,1);sum(Classes.H3,1);sum(Classes.H4,1);...
    sum(Classes.H5,1);sum(Classes.H6,1);sum(Classes.H7,1);sum(Classes.H8,1);sum(Classes.H9,1);sum(Classes.H10,1)];



%% Infection
Infec = [sum(Classes.I1,1);sum(Classes.I2,1);sum(Classes.I3,1);sum(Classes.I4,1);...
    sum(Classes.I5,1);sum(Classes.I6,1);sum(Classes.I7,1);sum(Classes.I8,1);sum(Classes.I9,1);sum(Classes.I10,1)];

%% Death

CaseF = [sum(De1,1);sum(De2,1);sum(De3,1);sum(De4,1);sum(De5,1);sum(De6,1);sum(De7,1);sum(De8,1);sum(De9,1);sum(De10,1)];

%% Total Infection

%% GROUP ONE

prev1 = sum((Classes.A1 + Classes.J1 + Classes.Q1 +Classes.I1 +Classes.H1),1);

%% GROUP TWO

prev2 = sum((Classes.A2 + Classes.J2 + Classes.Q2 +Classes.I2 +Classes.H2),1);

%% GROUP THREE

prev3 = sum((Classes.A3 + Classes.J3 + Classes.Q3 +Classes.I3 +Classes.H3),1);

%% GROUP FOUR

prev4 = sum((Classes.A4 + Classes.J4 + Classes.Q4 +Classes.I4 +Classes.H4),1);

%% GROUP FIVE

prev5 = sum((Classes.A5 + Classes.J5 + Classes.Q5 +Classes.I5 +Classes.H5),1);

%% GROUP SIX

prev6 = sum((Classes.A6 + Classes.J6 + Classes.Q6 +Classes.I6 +Classes.H6),1);

%% GROUP SEVEN

prev7 = sum((Classes.A7 + Classes.J7 + Classes.Q7 +Classes.I7 +Classes.H7),1);

%% GROUP EIGHT

prev8 = sum((Classes.A8 + Classes.J8 + Classes.Q8 +Classes.I8 +Classes.H8),1);

%% GROUP NINE

prev9 = sum((Classes.A9 + Classes.J9 + Classes.Q9 +Classes.I9 +Classes.H9),1);

%% GROUP TEN

prev10 = sum((Classes.A10 + Classes.J10 + Classes.Q10 +Classes.I10 +Classes.H10),1);

PREV2 = [prev1;prev2;prev3;prev4;prev5;prev6;prev7;prev8;prev9;prev10];

%% EXPORT THE DATA
% writematrix(Infec,'Agesocial.xlsx', 'Sheet', 'Infection')
% writematrix(PREV2,'Agesocial.xlsx', 'Sheet', 'DailyInfec')
% writematrix(Hospital,'Agesocial.xlsx', 'Sheet', 'Hospitalised')
% writematrix(CaseF,'Agesocial.xlsx', 'Sheet', 'Death')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% HEAT MAP COMBINED PLOT %%%%%%%%%%%%%%%%%%%%%

% clims = [min(Hospital(:)),max(Hospital(:))]; %Scale of the heatmap
customlabels2 = {'0-4', '5-9' ,'10-14', '15-19', '20-24', '25-29','30-34', ...
        '35-39','40-44','45-49','50-54','55-59','60-64','65-69','70-74','75-79', ...
        '80-84','85-89','90-94','95-99','99+'};
xValues = 1:21;
yValues = 1:10;

customlabels3 = {'1', '2' ,'3', '4', '5', '6','7','8','9','10'};
% figure;
% imagesc(PREV2,clims)
% colormap(flipud(hot))  %(hsv(512))
% colorbar
% title("Hospitalised")
% xticks(xValues);
% yticks(yValues);
% xticklabels(customlabels2);
% yticklabels(customlabels3);
% xlabel('Age Groups')
% ylabel('Deprivation Decile')
% grid on

figure;
subplot(2,2,1)
clims1 = [min(PREV2(:)),max(PREV2(:))]; %Scale of the heatmap
imagesc(PREV2,clims1)
colormap(flipud(hot))  %(hsv(512))
colorbar
title("(a) Daily infection.")
xticks(xValues);
yticks(yValues);
xticklabels(customlabels2);
yticklabels(customlabels3);
% xlabel('Age Groups')
ylabel('Deprivation Decile')
hold off

subplot(2,2,2)
clims2 = [min(Infec(:)),max(Infec(:))]; %Scale of the heatmap
imagesc(Infec,clims2)
colormap(flipud(hot))  %(hsv(512))
colorbar
title("(b) Reported daily infection")
xticks(xValues);
yticks(yValues);
xticklabels(customlabels2);
yticklabels(customlabels3);
% xlabel('Age Groups')
ylabel('Deprivation Decile')
hold off

subplot(2,2,3)
clims3 = [min(Hospital(:)),max(Hospital(:))]; %Scale of the heatmap
imagesc(Hospital,clims3)
colormap(flipud(hot))  %(hsv(512))
colorbar
title("(c) Hospitalised individuals")
xticks(xValues);
yticks(yValues);
xticklabels(customlabels2);
yticklabels(customlabels3);
xlabel('Age Groups')
ylabel('Deprivation Decile')
hold off

subplot(2,2,4)
clims4 = [min(CaseF(:)),max(CaseF(:))]; %Scale of the heatmap
imagesc(CaseF,clims4)
colormap(flipud(hot))  %(hsv(512))
colorbar
title("(d) Number of death")
xticks(xValues);
yticks(yValues);
xticklabels(customlabels2);
yticklabels(customlabels3);
xlabel('Age Groups')
ylabel('Deprivation Decile')
hold off






% figure;
% bar(Infec);
% title('Age Groups vs Social Groups');
% xlabel('Social Groups');
% ylabel('Values');


% surf(Infec);
% title('Social Groups vs Age Groups');
% xlabel('Age Groups');
% ylabel('Social Groups');
% zlabel('Values');

% clims = [min(Hospital(:)),max(Hospital(:))]; %Scale of the heatmap
% customlabels2 = {'0-4', '5-9' ,'10-14', '15-19', '20-24', '25-29','30-34', ...
%         '35-39','40-44','45-49','50-54','55-59','60-64','65-69','70-74','75-79', ...
%         '80-84','85-89','90-94','95-99','99+'};
% xValues = 1:21;
% yValues = 1:10;
% 
% customlabels3 = {'1', '2' ,'3', '4', '5', '6','7','8','9','10'};
% figure;
% imagesc(Hospital,clims)
% colormap(flipud(hot))  %(hsv(512))
% colorbar
% title("Hospitalised")
% xticks(xValues);
% yticks(yValues);
% xticklabels(customlabels2);
% yticklabels(customlabels3);
% xlabel('Age Groups')
% ylabel('Deprivation Decile')
% grid on

% 

% figure;
% imagesc(CaseF,[0,95000])
% colormap("turbo")  %(hsv(512))
% colorbar
% title("Case Fatality")
% xticks(xValues);
% yticks(yValues);
% xticklabels(customlabels2);
% yticklabels(customlabels3);
% xlabel('Age Groups')
% ylabel('Deprivation Decile')
% grid on
% 

% figure;
% imagesc(PREV2,[0,3800000])
% colormap("turbo")  %(hsv(512))
% colorbar
% title("TOTAL INFECTION")
% xticks(xValues);
% yticks(yValues);
% xticklabels(customlabels2);
% yticklabels(customlabels3);
% xlabel('Age Groups')
% ylabel('Deprivation Decile')
% grid on
% 
% 
% 
% %% Reported cases
% 
% IR2 = [sum(Classes.I1,1);sum(Classes.I2,1);sum(Classes.I3,1);sum(Classes.I4,1);...
%     sum(Classes.I5,1);sum(Classes.I6,1);sum(Classes.I7,1);sum(Classes.I8,1);...
%     sum(Classes.I9,1);sum(Classes.I10,1)];
% 
% figure;
% imagesc(IR2,[0,400000])
% colormap("turbo")  %(hsv(512))
% colorbar
% title("REPORTED CASES")
% xticks(xValues);
% yticks(yValues);
% xticklabels(customlabels2);
% yticklabels(customlabels3);
% xlabel('Age Groups')
% ylabel('Deprivation Decile')
% grid on
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %% COMBINE PLOT
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure;
% subplot(2,1,1)
% imagesc(PREV2,[min(PREV2(:)),max(PREV2(:))])
% colormap("turbo")  %(hsv(512))
% colorbar
% title("TOTAL INFECTION")
% set(gca,'xtick',[]);
% % xticks(xValues);
% yticks(yValues);
% % xticklabels(customlabels2);
% yticklabels(customlabels3);
% % xlabel('Age Groups')
% ylabel('Deprivation Decile')
% grid on
% 
% subplot(2,1,2)
% imagesc(IR2,[min(IR2(:)),max(IR2(:))])
% colormap("turbo")  %(hsv(512))
% colorbar
% title("REPORTED CASES")
% xticks(xValues);
% yticks(yValues);
% xticklabels(customlabels2);
% yticklabels(customlabels3);
% xlabel('Age Groups')
% ylabel('Deprivation Decile')
% grid on
% 
% figure;
% subplot(2,1,1)
% imagesc(Hospital,[min(Hospital(:)),max(Hospital(:))])
% colormap("turbo")  %(hsv(512))
% colorbar
% title("Hospitalised")
% set(gca,'xtick',[]);
% % xticks(xValues);
% yticks(yValues);
% % xticklabels(customlabels2);
% yticklabels(customlabels3);
% % xlabel('Age Groups')
% ylabel('Deprivation Decile')
% grid on
% 
% subplot(2,1,2)
% imagesc(CaseF,[min(CaseF(:)),max(CaseF(:))])
% colormap("turbo")  %(hsv(512))
% colorbar
% title("Case Fatality")
% xticks(xValues);
% yticks(yValues);
% xticklabels(customlabels2);
% yticklabels(customlabels3);
% xlabel('Age Groups')
% ylabel('Deprivation Decile')
% grid on
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 








% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % TotalKnown = Hosp + IR;
% % figure;
% % plot(Classes.t, TotalKnown/100000)
% % xlabel('Time (days)')
% % ylabel({' ``Known" Infections'; '(per 100,000 people)'})
% % legend("1",'2','3','4','5','6','7','8','9','10')
% % hold on 
% % rectangle('Position',[83 0 71 18/100],'FaceColor', [1, 0, 0, 0.2],'EdgeColor', [1, 0, 0, 0.2])
% % hold off
% % 
% % 
% % %% PLOT
% % 
% % figure(7)
% % plot(Classes.t,TOTAL_DE/100000)
% % xlabel('time (days)')
% % ylabel({'Daily Fatality Reported'; '(per 100,000 people)'})
% % legend("1",'2','3','4','5','6','7','8','9','10')
% % hold on 
% % rectangle('Position',[83 0 71 3500/100000],'FaceColor', [1, 0, 0, 0.2],'EdgeColor', [1, 0, 0, 0.2])
% % hold off
% % 
% % 
% % figure(8)
% % subplot(2,2,[1 3])
% % plot(Classes.t,Hosp/100000)
% % xlabel('time (days)')
% % ylabel({'Hospital Cases'; '(per 100,000 people)'})
% % legend("1",'2','3','4','5','6','7','8','9','10')
% % hold on 
% % rectangle('Position',[83 0 71 1600/10000],'FaceColor', [1, 0, 0, 0.2],'EdgeColor', [1, 0, 0, 0.2])
% % 
% % subplot(2,2,2)
% % plot(Classes.t,Death/100000)
% % xlabel('time (days)')
% % ylabel({'Fatally Infectious'; '(per 100,000 people)'})
% % legend("1",'2','3','4','5','6','7','8','9','10')
% % 
% % subplot(2,2,4)
% % plot(Classes.t,TOTAL_DE/100000)
% % xlabel('time (days)')
% % ylabel({'Daily Fatality Reported'; '(per 100,000 people)'})
% % legend("1",'2','3','4','5','6','7','8','9','10')
% % hold on 
% % rectangle('Position',[83 0 71 3000/100000],'FaceColor', [1, 0, 0, 0.2],'EdgeColor', [1, 0, 0, 0.2])
% % hold off
% % 

% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %% Hospitalised by age
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Hospital = {Classes.H1,Classes.H2,Classes.H3,Classes.H4,Classes.H5,...
% %     Classes.H6,Classes.H7,Classes.H8,Classes.H9,Classes.H10};
% % 
% % for i = 1:2:10
% %     % Create a figure for each pair of sets
% %     figure;
% % 
% %     % Plot the first set in the pair
% %     subplot(1, 2, 1);
% %     plot(Classes.t,Hospital{i}/100000); % You can adjust the layout, colormap, and other properties
% %     title(['Group ' num2str(i)]);
% %     ylabel({'Hospital Cases'; '(per 100,000 people)'})
% %     xlabel("time (days)")
% %     hold on
% %     rectangle('Position',[83 0 71 1/100],'FaceColor', [1, 0, 0, 0.2],'EdgeColor', [1, 0, 0, 0.2])
% % 
% %     % Plot the second set in the pair
% %     subplot(1, 2, 2);
% %     plot(Classes.t,Hospital{i}/100000); % You can adjust the layout, colormap, and other properties
% %     title(['Group ' num2str(i + 1)]);
% %     xlabel("time (days)")
% %     hold on
% %     legend('0-4', '5-9' ,'10-14', '15-19', '20-24', '25-29','30-34', ...
% %         '35-39','40-44','45-49','50-54','55-59','60-64','65-69','70-74','75-79', ...
% %         '80-84','85-89','90-94','95-99','99+')
% %     hold on
% %     rectangle('Position',[83 0 71 1/100],'FaceColor', [1, 0, 0, 0.2],'EdgeColor', [1, 0, 0, 0.2])
% % end
% % 
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %% NUMBER OF REPORTED CASES  age
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DailyInc = {Classes.I1,Classes.I2,Classes.I3,Classes.I4,Classes.I5,Classes.I6,Classes.I7,Classes.I8,Classes.I9,Classes.I10};
% for i = 1:2:10
%     % Create a figure for each pair of sets
%     figure;
% 
%     % Plot the first set in the pair
%     subplot(1, 2, 1);
%     plot(Classes.t,DailyInc{i}/100000); % You can adjust the layout, colormap, and other properties
%     title(['Group ' num2str(i)]);
%     ylabel({'Daily Reported Cases'; '(per 100,000 people)'})
%     xlabel("time (days)")
%     hold on
%     rectangle('Position',[83 0 71 3/1000],'FaceColor', [1, 0, 0, 0.2],'EdgeColor', [1, 0, 0, 0.2])
% 
%     % Plot the second set in the pair
%     subplot(1, 2, 2);
%     plot(Classes.t,DailyInc{i}/100000); % You can adjust the layout, colormap, and other properties
%     title(['Group ' num2str(i + 1)]);
%     xlabel("time (days)")
%     hold on
%     legend('0-4', '5-9' ,'10-14', '15-19', '20-24', '25-29','30-34', ...
%         '35-39','40-44','45-49','50-54','55-59','60-64','65-69','70-74','75-79', ...
%         '80-84','85-89','90-94','95-99','99+')
%     hold on
%     rectangle('Position',[83 0 71 3/1000],'FaceColor', [1, 0, 0, 0.2],'EdgeColor', [1, 0, 0, 0.2])
% end

% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%% Death BY SOCIODEMOGRAPHIC GROUP %%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % DEATH_AGE = {De1,De2,De3,De4,De5,De6,De7,De8,De9,De10};
% % 
% % for i = 1:2:10
% %     % Create a figure for each pair of sets
% %     figure;
% % 
% %     % Plot the first set in the pair
% %     subplot(1, 2, 1);
% %     plot(Classes.t,DEATH_AGE{i}/100000); % You can adjust the layout, colormap, and other properties
% %     title(['Group ' num2str(i)]);
% %     ylabel({'Case Fatality'; '(per 100,000 people)'})
% %     xlabel("time (days)")
% %     hold on
% %     rectangle('Position',[83 0 71 600/100000],'FaceColor', [1, 0, 0, 0.2],'EdgeColor', [1, 0, 0, 0.2])
% % 
% %     % Plot the second set in the pair
% %     subplot(1, 2, 2);
% %     plot(Classes.t,DEATH_AGE{i}/100000); % You can adjust the layout, colormap, and other properties
% %     title(['Group ' num2str(i + 1)]);
% %     xlabel("time (days)")
% %     hold on
% %     legend('0-4', '5-9' ,'10-14', '15-19', '20-24', '25-29','30-34', ...
% %         '35-39','40-44','45-49','50-54','55-59','60-64','65-69','70-74','75-79', ...
% %         '80-84','85-89','90-94','95-99','99+')
% %     hold on
% %     rectangle('Position',[83 0 71 600/100000],'FaceColor', [1, 0, 0, 0.2],'EdgeColor', [1, 0, 0, 0.2])
% % end
% % 
% % 
% % %% prevelance by age
% % Prev1 = (Classes.A1 + Classes.J1 + Classes.Q1 +Classes.I1 +Classes.H1);
% % Prev2 = (Classes.A2 + Classes.J2 + Classes.Q2 +Classes.I2 +Classes.H2);
% % Prev3 = (Classes.A3 + Classes.J3 + Classes.Q3 +Classes.I3 +Classes.H3);
% % Prev4 = (Classes.A4 + Classes.J4 + Classes.Q4 +Classes.I4 +Classes.H4);
% % Prev5 = (Classes.A5 + Classes.J5 + Classes.Q5 +Classes.I5 +Classes.H5);
% % Prev6 = (Classes.A6 + Classes.J6 + Classes.Q6 +Classes.I6 +Classes.H6);
% % Prev7 = (Classes.A7 + Classes.J7 + Classes.Q7 +Classes.I7 +Classes.H7);
% % Prev8 = (Classes.A8 + Classes.J8 + Classes.Q8 +Classes.I8 +Classes.H8);
% % Prev9 = (Classes.A9 + Classes.J9 + Classes.Q9 +Classes.I9 +Classes.H9);
% % Prev10 = (Classes.A10 + Classes.J10 + Classes.Q10 +Classes.I10 +Classes.H10);
% % 
% % PrevAge = {Prev1,Prev2,Prev3,Prev4,Prev5,Prev6,Prev7,Prev8,Prev9,Prev10};
% % 
% % for i = 1:2:10
% %     % Create a figure for each pair of sets
% %     figure;
% % 
% %     % Plot the first set in the pair
% %     subplot(1, 2, 1);
% %     plot(Classes.t,PrevAge{i}/100000); % You can adjust the layout, colormap, and other properties
% %     title(['Group ' num2str(i)]);
% %     ylabel({'Total Number Infected'; '(per 100,000 people)'})
% %     xlabel("time (days)")
% %     hold on
% %     rectangle('Position',[83 0 71 12/100],'FaceColor', [1, 0, 0, 0.2],'EdgeColor', [1, 0, 0, 0.2])
% % 
% %     % Plot the second set in the pair
% %     subplot(1, 2, 2);
% %     plot(Classes.t,PrevAge{i}/100000); % You can adjust the layout, colormap, and other properties
% %     title(['Group ' num2str(i + 1)]);
% %     xlabel("time (days)")
% %     hold on
% %     legend('0-4', '5-9' ,'10-14', '15-19', '20-24', '25-29','30-34', ...
% %         '35-39','40-44','45-49','50-54','55-59','60-64','65-69','70-74','75-79', ...
% %         '80-84','85-89','90-94','95-99','99+')
% %     hold on
% %     rectangle('Position',[83 0 71 12/100],'FaceColor', [1, 0, 0, 0.2],'EdgeColor', [1, 0, 0, 0.2])
% % end
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %% SUBSET TOTAL MAXMUM INFECTION TIME BY AGE
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % %%  PREVALENCE BY AGE
% % 
% % %% Before and during lockdown
% % start_index = find(Classes.t >= 0, 1);
% % end_index = find(Classes.t<= 153,1,"last");
% % 
% % %% prevelance by age
% % Prev1 = (Classes.A1 + Classes.J1 + Classes.Q1 +Classes.I1 +Classes.H1);
% % Age_prev1 = round(max(Prev1(start_index:end_index,:)));
% % 
% % Prev2 = (Classes.A2 + Classes.J2 + Classes.Q2 +Classes.I2 +Classes.H2);
% % Age_prev2 = round(max(Prev2(start_index:end_index,:)));
% % 
% % Prev3 = (Classes.A3 + Classes.J3 + Classes.Q3 +Classes.I3 +Classes.H3);
% % Age_prev3 = round(max(Prev3(start_index:end_index,:)));
% % 
% % Prev4 = (Classes.A4 + Classes.J4 + Classes.Q4 +Classes.I4 +Classes.H4);
% % Age_prev4 = round(max(Prev4(start_index:end_index,:)));
% % 
% % Prev5 = (Classes.A5 + Classes.J5 + Classes.Q5 +Classes.I5 +Classes.H5);
% % Age_prev5 = round(max(Prev5(start_index:end_index,:)));
% % 
% % Prev6 = (Classes.A6 + Classes.J6 + Classes.Q6 +Classes.I6 +Classes.H6);
% % Age_prev6 = round(max(Prev6(start_index:end_index,:)));
% % 
% % Prev7 = (Classes.A7 + Classes.J7 + Classes.Q7 +Classes.I7 +Classes.H7);
% % Age_prev7 = round(max(Prev7(start_index:end_index,:)));
% % 
% % Prev8 = (Classes.A8 + Classes.J8 + Classes.Q8 +Classes.I8 +Classes.H8);
% % Age_prev8 = round(max(Prev8(start_index:end_index,:)));
% % 
% % Prev9 = (Classes.A9 + Classes.J9 + Classes.Q9 +Classes.I9 +Classes.H9);
% % Age_prev9 = round(max(Prev9(start_index:end_index,:)));
% % 
% % Prev10 = (Classes.A10 + Classes.J10 + Classes.Q10 +Classes.I10 +Classes.H10);
% % Age_prev10 = round(max(Prev10(start_index:end_index,:)));
% % 
% % PrevAge = [Age_prev1;Age_prev2;Age_prev3;Age_prev4;Age_prev5;...
% %     Age_prev6;Age_prev7;Age_prev8;Age_prev9;Age_prev10];
% % 
% % 
% % figure(12)
% % hold on
% % for i = 1:10
% %     plot(PrevAge(i,:))
% %     set(gca,'xtick', 1:1:21)
% % end
% % hold off
% % legend("1",'2','3','4','5','6','7','8','9','10')
% % title("Before and during lockdown")
% % 
% % %% After Lockdown
% % start_index2 = find(Classes.t >= 154, 1);
% % end_index2 = find(Classes.t<= 350,1,"last");
% % 
% % %% prevelance by age
% % aPrev1 = (Classes.A1 + Classes.J1 + Classes.Q1 +Classes.I1 +Classes.H1);
% % AAge_prev1 = round(max(aPrev1(start_index2:end_index2,:)));
% % 
% % aPrev2 = (Classes.A2 + Classes.J2 + Classes.Q2 +Classes.I2 +Classes.H2);
% % AAge_prev2 = round(max(aPrev2(start_index2:end_index2,:)));
% % 
% % aPrev3 = (Classes.A3 + Classes.J3 + Classes.Q3 +Classes.I3 +Classes.H3);
% % AAge_prev3 = round(max(aPrev3(start_index2:end_index2,:)));
% % 
% % aPrev4 = (Classes.A4 + Classes.J4 + Classes.Q4 +Classes.I4 +Classes.H4);
% % AAge_prev4 = round(max(aPrev4(start_index2:end_index2,:)));
% % 
% % aPrev5 = (Classes.A5 + Classes.J5 + Classes.Q5 +Classes.I5 +Classes.H5);
% % AAge_prev5 = round(max(aPrev5(start_index2:end_index2,:)));
% % 
% % aPrev6 = (Classes.A6 + Classes.J6 + Classes.Q6 +Classes.I6 +Classes.H6);
% % AAge_prev6 = round(max(aPrev6(start_index2:end_index2,:)));
% % 
% % aPrev7 = (Classes.A7 + Classes.J7 + Classes.Q7 +Classes.I7 +Classes.H7);
% % AAge_prev7 = round(max(aPrev7(start_index2:end_index2,:)));
% % 
% % aPrev8 = (Classes.A8 + Classes.J8 + Classes.Q8 +Classes.I8 +Classes.H8);
% % AAge_prev8 = round(max(aPrev8(start_index2:end_index2,:)));
% % 
% % aPrev9 = (Classes.A9 + Classes.J9 + Classes.Q9 +Classes.I9 +Classes.H9);
% % AAge_prev9 = round(max(aPrev9(start_index2:end_index2,:)));
% % 
% % aPrev10 = (Classes.A10 + Classes.J10 + Classes.Q10 +Classes.I10 +Classes.H10);
% % AAge_prev10 = round(max(aPrev10(start_index2:end_index2,:)));
% % 
% % PrevAge2 = [AAge_prev1;AAge_prev2;AAge_prev3;AAge_prev4;AAge_prev5;...
% %     AAge_prev6;AAge_prev7;AAge_prev8;AAge_prev9;AAge_prev10];
% % 
% % figure(13)
% % hold on
% % for i = 1:10
% %     plot(PrevAge2(i,:))
% %     set(gca,'xtick', 1:1:21)
% % end
% % hold off
% % legend("1",'2','3','4','5','6','7','8','9','10')
% % title("After lockdown")
% % 
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %% DEATH BY AGE GROUP
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % % Before and lockdown
% % 
% % Age_death1 = round(max(Classes.D1(start_index:end_index,:)));
% % Age_death2 = round(max(Classes.D2(start_index:end_index,:)));
% % Age_death3 = round(max(Classes.D3(start_index:end_index,:)));
% % Age_death4 = round(max(Classes.D4(start_index:end_index,:)));
% % Age_death5 = round(max(Classes.D5(start_index:end_index,:)));
% % Age_death6 = round(max(Classes.D6(start_index:end_index,:)));
% % Age_death7 = round(max(Classes.D7(start_index:end_index,:)));
% % Age_death8 = round(max(Classes.D8(start_index:end_index,:)));
% % Age_death9 = round(max(Classes.D9(start_index:end_index,:)));
% % Age_death10 = round(max(Classes.D10(start_index:end_index,:)));
% % 
% % Age_Death = [Age_death1;Age_death2;Age_death3;Age_death4;Age_death5;...
% %     Age_death6;Age_death7;Age_death8;Age_death9;Age_death10];
% % % After lockdown
% % 
% % Age_dea1 = round(max(Classes.D1(start_index2:end_index2,:)));
% % Age_dea2 = round(max(Classes.D2(start_index2:end_index2,:)));
% % Age_dea3 = round(max(Classes.D3(start_index2:end_index2,:)));
% % Age_dea4 = round(max(Classes.D4(start_index2:end_index2,:)));
% % Age_dea5 = round(max(Classes.D5(start_index2:end_index2,:)));
% % Age_dea6 = round(max(Classes.D6(start_index2:end_index2,:)));
% % Age_dea7 = round(max(Classes.D7(start_index2:end_index2,:)));
% % Age_dea8 = round(max(Classes.D8(start_index2:end_index2,:)));
% % Age_dea9 = round(max(Classes.D9(start_index2:end_index2,:)));
% % Age_dea10 = round(max(Classes.D10(start_index2:end_index2,:)));
% % 
% % Age_DeaD = [Age_dea1;Age_dea2;Age_dea3;Age_dea4;Age_dea5;...
% %     Age_dea6;Age_dea7;Age_dea8;Age_dea9;Age_dea10];
% % 
% % figure(14)
% % hold on
% % for i = 1:10
% %     plot(Age_Death(i,:))
% %     set(gca,'xtick', 1:1:21)
% % end
% % hold off
% % legend("1",'2','3','4','5','6','7','8','9','10')
% % title("Before and during lockdown")
% % 
% % 
% % figure(15)
% % hold on
% % for i = 1:10
% %     plot(Age_DeaD(i,:))
% %     set(gca,'xtick', 1:1:21)
% % end
% % hold off
% % legend("1",'2','3','4','5','6','7','8','9','10')
% % title("After lockdown")
% % 
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %% HOSPITALISATION BY AGE GROUP
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % %% Before and lockdown
% % 
% % Age_Hos1 = round(max(Classes.H1(start_index:end_index,:)));
% % Age_Hos2 = round(max(Classes.H2(start_index:end_index,:)));
% % Age_Hos3 = round(max(Classes.H3(start_index:end_index,:)));
% % Age_Hos4 = round(max(Classes.H4(start_index:end_index,:)));
% % Age_Hos5 = round(max(Classes.H5(start_index:end_index,:)));
% % Age_Hos6 = round(max(Classes.H6(start_index:end_index,:)));
% % Age_Hos7 = round(max(Classes.H7(start_index:end_index,:)));
% % Age_Hos8 = round(max(Classes.H8(start_index:end_index,:)));
% % Age_Hos9 = round(max(Classes.H9(start_index:end_index,:)));
% % Age_Hos10 = round(max(Classes.H10(start_index:end_index,:)));
% % 
% % Age_HOS = [Age_Hos1;Age_Hos2;Age_Hos3;Age_Hos4;Age_Hos5;...
% %     Age_Hos6;Age_Hos7;Age_Hos8;Age_Hos9;Age_Hos10];
% % 
% % %% After lockdown
% % 
% % Age_Hosp1 = round(max(Classes.H1(start_index2:end_index2,:)));
% % Age_Hosp2 = round(max(Classes.H2(start_index2:end_index2,:)));
% % Age_Hosp3 = round(max(Classes.H3(start_index2:end_index2,:)));
% % Age_Hosp4 = round(max(Classes.H4(start_index2:end_index2,:)));
% % Age_Hosp5 = round(max(Classes.H5(start_index2:end_index2,:)));
% % Age_Hosp6 = round(max(Classes.H6(start_index2:end_index2,:)));
% % Age_Hosp7 = round(max(Classes.H7(start_index2:end_index2,:)));
% % Age_Hosp8 = round(max(Classes.H8(start_index2:end_index2,:)));
% % Age_Hosp9 = round(max(Classes.H9(start_index2:end_index2,:)));
% % Age_Hosp10 = round(max(Classes.H10(start_index2:end_index2,:)));
% % 
% % Age_HOSP = [Age_Hosp1;Age_Hosp2;Age_Hosp3;Age_Hosp4;Age_Hosp5;...
% %     Age_Hosp6;Age_Hosp7;Age_Hosp8;Age_Hosp9;Age_Hosp10];
% % 
% % figure(16)
% % hold on
% % for i = 1:10
% %     plot(Age_HOS(i,:))
% %     set(gca,'xtick', 1:1:21)
% % end
% % hold off
% % legend("1",'2','3','4','5','6','7','8','9','10')
% % title("Before and during lockdown")
% % 
% % 
% % figure(17)
% % hold on
% % for i = 1:10
% %     plot(Age_HOSP(i,:))
% %     set(gca,'xtick', 1:1:21)
% % end
% % hold off
% % legend("1",'2','3','4','5','6','7','8','9','10')
% % title("After lockdown")
% % 
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %% DAILY INCIDENCE BY AGE GROUP
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % %% Before and lockdown
% % 
% % Age_inc1 = round(max(Inc1(start_index:end_index,:)));
% % Age_inc2 = round(max(Inc2(start_index:end_index,:)));
% % Age_inc3 = round(max(Inc3(start_index:end_index,:)));
% % Age_inc4 = round(max(Inc4(start_index:end_index,:)));
% % Age_inc5 = round(max(Inc5(start_index:end_index,:)));
% % Age_inc6 = round(max(Inc6(start_index:end_index,:)));
% % Age_inc7 = round(max(Inc7(start_index:end_index,:)));
% % Age_inc8 = round(max(Inc8(start_index:end_index,:)));
% % Age_inc9 = round(max(Inc9(start_index:end_index,:)));
% % Age_inc10 = round(max(Inc10(start_index:end_index,:)));
% % 
% % Age_INC = [Age_inc1;Age_inc2;Age_inc3;Age_inc4;Age_inc5;...
% %     Age_inc6;Age_inc7;Age_inc8;Age_inc9;Age_inc10];
% % 
% % %% After lockdown
% % 
% % Age_INCD1 = round(max(Inc1(start_index2:end_index2,:)));
% % Age_INCD2 = round(max(Inc2(start_index2:end_index2,:)));
% % Age_INCD3 = round(max(Inc3(start_index2:end_index2,:)));
% % Age_INCD4 = round(max(Inc4(start_index2:end_index2,:)));
% % Age_INCD5 = round(max(Inc5(start_index2:end_index2,:)));
% % Age_INCD6 = round(max(Inc6(start_index2:end_index2,:)));
% % Age_INCD7 = round(max(Inc7(start_index2:end_index2,:)));
% % Age_INCD8 = round(max(Inc8(start_index2:end_index2,:)));
% % Age_INCD9 = round(max(Inc9(start_index2:end_index2,:)));
% % Age_INCD10 = round(max(Inc10(start_index2:end_index2,:)));
% % 
% % Age_INCD = [Age_INCD1;Age_INCD2;Age_INCD3;Age_INCD4;Age_INCD5;...
% %     Age_INCD6;Age_INCD7;Age_INCD8;Age_INCD9;Age_INCD10];
% % 
% % figure(18)
% % hold on
% % for i = 1:10
% %     plot(Age_INC(i,:))
% %     set(gca,'xtick', 1:1:21)
% % end
% % hold off
% % legend("1",'2','3','4','5','6','7','8','9','10')
% % title("Before and during lockdown")
% % 
% % 
% % figure(19)
% % hold on
% % for i = 1:10
% %     plot(Age_INCD(i,:))
% %     set(gca,'xtick', 1:1:21)
% % end
% % hold off
% % legend("1",'2','3','4','5','6','7','8','9','10')
% % title("After lockdown")
% % 
% % 
% 
% 
% 
% 
% % figure(20)
% % for i = 1:10
% %     % Create subplots in a 2x5 grid 
% %     subplot(2, 5, i);
% % 
% %     % Plot each set of data
% %     plot(Classes.t,Hospital{i});
% %     title("Group",i)
% %     xlabel('time (days)')
% %     hold on
% % end
% % legend('0-4', '5-9' ,'10-14', '15-19', '20-24', '25-29','30-34', ...
% %         '35-39','40-44','45-49','50-54','55-59','60-64','65-69','70-74','75-79', ...
% %         '80-84','85-89','90-94','95-99','99+')
% % 
% % figure(21)
% % plot(Classes.t,Classes.H1)
% % xlabel('time (days)')
% % ylabel({'Hospital Cases'; '(per 100,000 people)'})
% % legend("1",'2','3','4','5','6','7','8','9','10')
% % hold off
% % 
% % 
% % 
% 
% 
% 
% 
% 
% 
% 
% % plot(Classes.t,prev1/100000)
% % hold on
% % plot(Classes.t,prev2/100000)
% % hold on
% % plot(Classes.t,prev3/100000)
% % hold on
% % plot(Classes.t,prev4/100000)
% % hold on
% % plot(Classes.t,prev5/100000)
% % hold on
% % plot(Classes.t,prev6/100000)
% % hold on
% % plot(Classes.t,prev7/100000)
% % hold on
% % plot(Classes.t,prev8/100000)
% % hold on
% % plot(Classes.t,prev9/100000)
% % hold on
% % plot(Classes.t,prev10/100000)
% % hold on
% 
% 
% 
% 
% 
% % hold on
% % plot(Classes.t,Total_Inc2/100000)
% % hold on
% % plot(Classes.t,Total_Inc3/100000)
% % hold on
% % plot(Classes.t,Total_Inc4/100000)
% % hold on
% % plot(Classes.t,Total_Inc5/100000)
% % hold on
% % plot(Classes.t,Total_Inc6/100000)
% % hold on
% % plot(Classes.t,Total_Inc7/100000)
% % hold on
% % plot(Classes.t,Total_Inc8/100000)
% % hold on
% % plot(Classes.t,Total_Inc9/100000)
% % hold on
% % plot(Classes.t,Total_Inc10/100000)
% % hold on
% 
% 
% 
% 
% % hold on
% % plot(Classes.t,sum(Classes.Ir2,2)/100000)
% % hold on
% % plot(Classes.t,sum(Classes.Ir3,2)/100000)
% % hold on
% % plot(Classes.t,sum(Classes.Ir4,2)/100000)
% % hold on
% % plot(Classes.t,sum(Classes.Ir5,2)/100000)
% % hold on
% % plot(Classes.t,sum(Classes.Ir6,2)/100000)
% % hold on
% % plot(Classes.t,sum(Classes.Ir7,2)/100000)
% % hold on
% % plot(Classes.t,sum(Classes.Ir8,2)/100000)
% % hold on
% % plot(Classes.t,sum(Classes.Ir9,2)/100000)
% % hold on
% % plot(Classes.t,sum(Classes.Ir10,2)/100000)
% % hold on
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%% NUMBER OF PEOPLE REPORTED PER DAY BY SOCIODEMOGRAPHIC GROUP %%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % %% Group One
% % 
% % Inc1 = [zeros(1,n);diff(Classes.Ir1)];
% % Total_Inc1 = sum(Inc1,2);
% % 
% % %% Group Two
% % 
% % Inc2 = [zeros(1,n);diff(Classes.Ir2)];
% % Total_Inc2 = sum(Inc2,2);
% % 
% % %% Group Three
% % 
% % Inc3 = [zeros(1,n);diff(Classes.Ir3)];
% % Total_Inc3 = sum(Inc3,2);
% % 
% % %% Group Four
% % 
% % Inc4 = [zeros(1,n);diff(Classes.Ir4)];
% % Total_Inc4 = sum(Inc4,2);
% % 
% % %% Group Five
% % 
% % Inc5 = [zeros(1,n);diff(Classes.Ir5)];
% % Total_Inc5 = sum(Inc5,2);
% % 
% % %% Group Six
% % 
% % Inc6 = [zeros(1,n);diff(Classes.Ir6)];
% % Total_Inc6 = sum(Inc6,2);
% % 
% % %% Group Seven
% % 
% % Inc7 = [zeros(1,n);diff(Classes.Ir7)];
% % Total_Inc7 = sum(Inc7,2);
% % 
% % %% Group Eight
% % 
% % Inc8 = [zeros(1,n);diff(Classes.Ir8)];
% % Total_Inc8 = sum(Inc8,2);
% % 
% % %% Group Nine
% % 
% % Inc9 = [zeros(1,n);diff(Classes.Ir9)];
% % Total_Inc9 = sum(Inc9,2);
% % 
% % %% Group Ten
% % 
% % Inc10 = [zeros(1,n);diff(Classes.Ir10)];
% % Total_Inc10 = sum(Inc10,2);
% % 
% % TOTAL_INC = [Total_Inc1,Total_Inc2,Total_Inc3,Total_Inc4,Total_Inc5,...
% %     Total_Inc6,Total_Inc7,Total_Inc8,Total_Inc9,Total_Inc10];
% % 
% % 
% % %% PLOT
% % 
% % figure(2)
% % plot(Classes.t,TOTAL_INC/100000)
% % xlabel('time (days)')
% % ylabel({'Daily Number of Reported Cases'; '(per 100,000 people)'})
% % legend("1",'2','3','4','5','6','7','8','9','10')
% % hold on 
% % rectangle('Position',[83 0 71 3500/100000],'FaceColor', [1, 0, 0, 0.2],'EdgeColor', [1, 0, 0, 0.2])
% % hold off


% 
% % figure;
% % subplot(2,1,1)
% % plot(Classes.t,PREV/100000)
% % xlabel('Time (days)')
% % ylabel({'Total Number of Infectious'; '(per 100,000 people)'})
% % legend("1",'2','3','4','5','6','7','8','9','10')
% % hold on 
% % rectangle('Position',[83 0 71 100000/100000],'FaceColor', [1, 0, 0, 0.2],'EdgeColor', [1, 0, 0, 0.2])
% % hold off
% % 
% % subplot(2,1,2)
% % plot(Classes.t,IR/100000)
% % xlabel('Time (days)')
% % ylabel({'Number of Reported Cases'; '(per 100,000 people)'})
% % legend("1",'2','3','4','5','6','7','8','9','10')
% % hold on 
% % rectangle('Position',[83 0 71 1/10],'FaceColor', [1, 0, 0, 0.2],'EdgeColor', [1, 0, 0, 0.2])
% % hold off
% % 
