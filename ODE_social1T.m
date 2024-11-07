%% The code for my PhD Project, this contains the age and social mixing
% Using a ODE15s for this analysis. 

%% Input 
% para - structure of paraeter
%ICs - Structure of initial conditions

%Mix_H - Age ixing in the house
%Mix_W - Age ixing in at work
%Mix_S - Age ixing in school
%Mix_O - Age ixing in other places
%beta - Total ixing that occur in all the locations
%N - Population size
% n- number of age groups
% soc - Soc demographic groups mixing matrix
% x - vector of mixing parttern


%% Output

%Classes - array of number of individual in each compartment

%% FUNCTION DEFINITION

function [Classes,R0,R02] = ODE_social1T(para,ICs,maxtime,Mix_H,Mix_W,Mix_S,Mix_O,n,soc,x)
%Set tolerance
opts = odeset('RelTol',1e-4,'AbsTol',1e-3); %tolerance level%tolerance level
     tspan = [0:1:maxtime]; %time span


%Shape of the initial condition
 All = reshape([ICs.S1'; ICs.E1'; ICs.A1'; ICs.J1';ICs.Q1';ICs.I1'; ICs.H1'; ICs.R1'; ICs.D1'; ICs.Ir1';...
                 ICs.S2'; ICs.E2'; ICs.A2'; ICs.J2';ICs.Q2';ICs.I2';ICs.H2'; ICs.R2'; ICs.D2'; ICs.Ir2'; ...
                 ICs.S3'; ICs.E3'; ICs.A3'; ICs.J3';ICs.Q3';ICs.I3';ICs.H3'; ICs.R3'; ICs.D3'; ICs.Ir3'; ...
                 ICs.S4'; ICs.E4'; ICs.A4'; ICs.J4';ICs.Q4';ICs.I4';ICs.H4'; ICs.R4'; ICs.D4'; ICs.Ir4'; ...
                 ICs.S5'; ICs.E5'; ICs.A5'; ICs.J5';ICs.Q5';ICs.I5';ICs.H5'; ICs.R5'; ICs.D5'; ICs.Ir5'; ...
                 ICs.S6'; ICs.E6'; ICs.A6'; ICs.J6';ICs.Q6';ICs.I6';ICs.H6'; ICs.R6'; ICs.D6'; ICs.Ir6'; ...
                 ICs.S7'; ICs.E7'; ICs.A7'; ICs.J7';ICs.Q7';ICs.I7';ICs.H7'; ICs.R7'; ICs.D7'; ICs.Ir7'; ...
                 ICs.S8'; ICs.E8'; ICs.A8'; ICs.J8';ICs.Q8';ICs.I8';ICs.H8'; ICs.R8'; ICs.D8'; ICs.Ir8'; ...
                 ICs.S9'; ICs.E9'; ICs.A9'; ICs.J9';ICs.Q9';ICs.I9';ICs.H9'; ICs.R9'; ICs.D9'; ICs.Ir9'; ...
                 ICs.S10'; ICs.E10'; ICs.A10'; ICs.J10';ICs.Q10';ICs.I10';ICs.H10';ICs.R10'; ICs.D10'; ICs.Ir10'], []*n,1);

    [t , pop] = ode15s(@(t,y)diff_socialmodel(t,y,para),tspan,All,opts);
     

 Classes = struct('S1',pop(:,1:1:n),'E1',pop(:,n+1:1:2*n),'A1',pop(:,(2*n)+1:1:3*n),'J1',pop(:,(3*n)+1:1:4*n),'Q1',pop(:,(4*n)+1:1:5*n),'I1',pop(:,(5*n)+1:1:6*n), ...
          'H1',pop(:,(6*n)+1:1:7*n),'R1',pop(:,(7*n)+1:1:8*n),'D1',pop(:,(8*n)+1:1:9*n),'Ir1',pop(:,(9*n)+1:1:10*n), ...
          'S2',pop(:,(10*n)+1:1:11*n),'E2',pop(:,(11*n)+1:1:12*n),'A2',pop(:,(12*n)+1:1:13*n),'J2',pop(:,(13*n)+1:1:14*n),'Q2',pop(:,(14*n)+1:1:15*n), ...
          'I2',pop(:,(15*n)+1:1:16*n),'H2',pop(:,(16*n)+1:1:17*n),'R2',pop(:,(17*n)+1:1:18*n),'D2', pop(:,(18*n)+1:1:19*n),'Ir2', pop(:,(19*n)+1:1:20*n), ...
          'S3',pop(:,(20*n)+1:1:21*n),'E3',pop(:,(21*n)+1:1:22*n),'A3',pop(:,(22*n)+1:1:23*n),'J3',pop(:,(23*n)+1:1:24*n),'Q3',pop(:,(24*n)+1:1:25*n), ...
          'I3',pop(:,(25*n)+1:1:26*n),'H3',pop(:,(26*n)+1:1:27*n),'R3',pop(:,(27*n)+1:1:28*n),'D3', pop(:,(28*n)+1:1:29*n),'Ir3', pop(:,(29*n)+1:1:30*n), ...
          'S4',pop(:,(30*n)+1:1:31*n),'E4',pop(:,(31*n)+1:1:32*n),'A4',pop(:,(32*n)+1:1:33*n),'J4',pop(:,(33*n)+1:1:34*n),'Q4',pop(:,(34*n)+1:1:35*n), ...
          'I4',pop(:,(35*n)+1:1:36*n),'H4',pop(:,(36*n)+1:1:37*n),'R4',pop(:,(37*n)+1:1:38*n),'D4', pop(:,(38*n)+1:1:39*n),'Ir4',pop(:,(39*n)+1:1:40*n), ...
          'S5',pop(:,(40*n)+1:1:41*n),'E5',pop(:,(41*n)+1:1:42*n),'A5',pop(:,(42*n)+1:1:43*n),'J5',pop(:,(43*n)+1:1:44*n),'Q5',pop(:,(44*n)+1:1:45*n), ...
          'I5',pop(:,(45*n)+1:1:46*n),'H5',pop(:,(46*n)+1:1:47*n),'R5',pop(:,(47*n)+1:1:48*n),'D5', pop(:,(48*n)+1:1:49*n),'Ir5', pop(:,(49*n)+1:1:50*n), ...
          'S6',pop(:,(50*n)+1:1:51*n),'E6',pop(:,(51*n)+1:1:52*n),'A6',pop(:,(52*n)+1:1:53*n),'J6',pop(:,(53*n)+1:1:54*n),'Q6',pop(:,(54*n)+1:1:55*n), ...
          'I6',pop(:,(55*n)+1:1:56*n),'H6',pop(:,(56*n)+1:1:57*n),'R6',pop(:,(57*n)+1:1:58*n),'D6', pop(:,(58*n)+1:1:59*n),'Ir6', pop(:,(59*n)+1:1:60*n), ...
          'S7',pop(:,(60*n)+1:1:61*n),'E7',pop(:,(61*n)+1:1:62*n),'A7',pop(:,(62*n)+1:1:63*n),'J7',pop(:,(63*n)+1:1:64*n),'Q7',pop(:,(64*n)+1:1:65*n), ...
          'I7',pop(:,(65*n)+1:1:66*n),'H7',pop(:,(66*n)+1:1:67*n),'R7',pop(:,(67*n)+1:1:68*n),'D7', pop(:,(68*n)+1:1:69*n),'Ir7', pop(:,(69*n)+1:1:70*n), ...
          'S8',pop(:,(70*n)+1:1:71*n),'E8',pop(:,(71*n)+1:1:72*n),'A8',pop(:,(72*n)+1:1:73*n),'J8',pop(:,(73*n)+1:1:74*n),'Q8',pop(:,(74*n)+1:1:75*n), ...
          'I8',pop(:,(75*n)+1:1:76*n),'H8',pop(:,(76*n)+1:1:77*n),'R8',pop(:,(77*n)+1:1:78*n),'D8', pop(:,(78*n)+1:1:79*n),'Ir8', pop(:,(79*n)+1:1:80*n), ...
          'S9',pop(:,(80*n)+1:1:81*n),'E9',pop(:,(81*n)+1:1:82*n),'A9',pop(:,(82*n)+1:1:83*n),'J9',pop(:,(83*n)+1:1:84*n),'Q9',pop(:,(84*n)+1:1:85*n), ...
          'I9',pop(:,(85*n)+1:1:86*n),'H9',pop(:,(86*n)+1:1:87*n),'R9',pop(:,(87*n)+1:1:88*n),'D9', pop(:,(88*n)+1:1:89*n),'Ir9', pop(:,(89*n)+1:1:90*n), ...
          'S10',pop(:,(90*n)+1:1:91*n),'E10',pop(:,(91*n)+1:1:92*n),'A10',pop(:,(92*n)+1:1:93*n),'J10',pop(:,(93*n)+1:1:94*n),'Q10',pop(:,(94*n)+1:1:95*n), ...
          'I10',pop(:,(95*n)+1:1:96*n),'H10',pop(:,(96*n)+1:1:97*n),'R10',pop(:,(97*n)+1:1:98*n),'D10', pop(:,(98*n)+1:1:99*n),'Ir10', pop(:,(99*n)+1:1:100*n),'t',t);

  %% Set up the population change

 function dpop = diff_socialmodel(t,pop,para)

      %% Run the logistic curve
      y = GeneralisedRichard(maxtime);
                    
                    for m =1:length(y)
                        [pi(m)] = y(m)/sum(sum(para.N));
                        
                    end

           
            %% Estimate the R0 value
             AgeMat0 =(Mix_H + Mix_W + Mix_S + Mix_O)*(para.a.*para.h'); % Define the mixing matrix at time t=0
             AS_0 = kron(soc,AgeMat0);
     % CREATE A CELL ARRAY FOR EACH ROW

     AgeMSoc = {[AS_0(1:21,1:21)],[AS_0(1:21,22:42)],[AS_0(1:21,43:63)],[AS_0(1:21,64:84)],[AS_0(1:21,85:105)],[AS_0(1:21,106:126)],[AS_0(1:21,127:147)],[AS_0(1:21,148:168)],[AS_0(1:21,169:189)],[AS_0(1:21,190:210)];
               [AS_0(22:42,1:21)],[AS_0(22:42,22:42)],[AS_0(22:42,43:63)],[AS_0(22:42,64:84)],[AS_0(22:42,85:105)],[AS_0(22:42,106:126)],[AS_0(22:42,127:147)],[AS_0(22:42,148:168)],[AS_0(22:42,169:189)],[AS_0(22:42,190:210)];
               [AS_0(43:63,1:21)],[AS_0(43:63,22:42)],[AS_0(43:63,43:63)],[AS_0(43:63,64:84)],[AS_0(43:63,85:105)],[AS_0(43:63,106:126)],[AS_0(43:63,127:147)],[AS_0(43:63,148:168)],[AS_0(43:63,169:189)],[AS_0(43:63,190:210)];
               [AS_0(64:84,1:21)],[AS_0(64:84,22:42)],[AS_0(64:84,43:63)],[AS_0(64:84,64:84)],[AS_0(64:84,85:105)],[AS_0(64:84,106:126)],[AS_0(64:84,127:147)],[AS_0(64:84,148:168)],[AS_0(64:84,169:189)],[AS_0(64:84,190:210)];
               [AS_0(85:105,1:21)],[AS_0(85:105,22:42)],[AS_0(85:105,43:63)],[AS_0(85:105,64:84)],[AS_0(85:105,85:105)],[AS_0(85:105,106:126)],[AS_0(85:105,127:147)],[AS_0(85:105,148:168)],[AS_0(85:105,169:189)],[AS_0(85:105,190:210)];
               [AS_0(106:126,1:21)],[AS_0(106:126,22:42)],[AS_0(106:126,43:63)],[AS_0(106:126,64:84)],[AS_0(106:126,85:105)],[AS_0(106:126,106:126)],[AS_0(106:126,127:147)],[AS_0(106:126,148:168)],[AS_0(106:126,169:189)],[AS_0(106:126,190:210)];
               [AS_0(127:147,1:21)],[AS_0(127:147,22:42)],[AS_0(127:147,43:63)],[AS_0(127:147,64:84)],[AS_0(127:147,85:105)],[AS_0(127:147,106:126)],[AS_0(127:147,127:147)],[AS_0(127:147,148:168)],[AS_0(127:147,169:189)],[AS_0(127:147,190:210)];
               [AS_0(148:168,1:21)],[AS_0(148:168,22:42)],[AS_0(148:168,43:63)],[AS_0(148:168,64:84)],[AS_0(148:168,85:105)],[AS_0(148:168,106:126)],[AS_0(148:168,127:147)],[AS_0(148:168,148:168)],[AS_0(148:168,169:189)],[AS_0(148:168,190:210)];
               [AS_0(169:189,1:21)],[AS_0(169:189,22:42)],[AS_0(169:189,43:63)],[AS_0(169:189,64:84)],[AS_0(169:189,85:105)],[AS_0(169:189,106:126)],[AS_0(169:189,127:147)],[AS_0(169:189,148:168)],[AS_0(169:189,169:189)],[AS_0(169:189,190:210)];
               [AS_0(190:210,1:21)],[AS_0(190:210,22:42)],[AS_0(190:210,43:63)],[AS_0(190:210,64:84)],[AS_0(190:210,85:105)],[AS_0(190:210,106:126)],[AS_0(190:210,127:147)],[AS_0(190:210,148:168)],[AS_0(190:210,169:189)],[AS_0(190:210,190:210)]};

      % Preallocate a cell array to store the maximum eigenvalues
                    maxeig = cell(size( AgeMSoc ));
                    R0 = cell(size( AgeMSoc ));

                    % Iterate over each element in AgeSoc and calculate the maximum eigenvalue

                    for i = 1:size( AgeMSoc , 1)
                        for j = 1:size( AgeMSoc , 2)
                            

                             % AgeMSoc {i, j} =  AgeMSoc {i, j} ;
                             
                            maxeig{i, j} = max(eig(  AgeMSoc{i, j}));
                            riskgroups(i,j) = maxeig{i, j};
                            
                 
                        end
                    end

                    R0= (((riskgroups* para.beta) / para.gamma)*(1-pi(1))*para.rho);
                    
                    R02= max(eig(R0));


 
             %% Capture the impact of different COVID-19 strain on infectivity
        
               
        
                         if (t > 0) && (t <= 83)
                                ageMix = Mix_H+Mix_S+Mix_W+Mix_O;
                                w = soc;
                           elseif (t > 83) && (t <= 164)
                               ageMix = NPIchange(Mix_H,Mix_S,Mix_W,Mix_O,"TOTAL_LOCKDOWN");
                               w = 0.5*(soc); 
                           elseif (t > 164) && (t <= 229)
                               ageMix = NPIchange(Mix_H,Mix_S,Mix_W,Mix_O,"EASING");
                                w = (soc); 
                            elseif (t > 229) && (t <= 290)
                               ageMix = NPIchange(Mix_H,Mix_S,Mix_W,Mix_O,"RELAXED_RESTRICTION");
                                w = (soc); 
                            else    
                               ageMix = NPIchange(Mix_H,Mix_S,Mix_W,Mix_O,"RELAXED_RESTRICTION");
                                w = (soc);
                         end



 % ageMix = Mix_H+Mix_S+Mix_W+Mix_O;
 % w = soc;
            
     %% Define Age Mixing

    
     AgeMat = ageMix*(para.a.*para.h');
     AS = kron(w,AgeMat);
     % CREATE A CELL ARRAY FOR EACH ROW

     AgeSoc = {[AS(1:21,1:21)],[AS(1:21,22:42)],[AS(1:21,43:63)],[AS(1:21,64:84)],[AS(1:21,85:105)],[AS(1:21,106:126)],[AS(1:21,127:147)],[AS(1:21,148:168)],[AS(1:21,169:189)],[AS(1:21,190:210)];
               [AS(22:42,1:21)],[AS(22:42,22:42)],[AS(22:42,43:63)],[AS(22:42,64:84)],[AS(22:42,85:105)],[AS(22:42,106:126)],[AS(22:42,127:147)],[AS(22:42,148:168)],[AS(22:42,169:189)],[AS(22:42,190:210)];
               [AS(43:63,1:21)],[AS(43:63,22:42)],[AS(43:63,43:63)],[AS(43:63,64:84)],[AS(43:63,85:105)],[AS(43:63,106:126)],[AS(43:63,127:147)],[AS(43:63,148:168)],[AS(43:63,169:189)],[AS(43:63,190:210)];
               [AS(64:84,1:21)],[AS(64:84,22:42)],[AS(64:84,43:63)],[AS(64:84,64:84)],[AS(64:84,85:105)],[AS(64:84,106:126)],[AS(64:84,127:147)],[AS(64:84,148:168)],[AS(64:84,169:189)],[AS(64:84,190:210)];
               [AS(85:105,1:21)],[AS(85:105,22:42)],[AS(85:105,43:63)],[AS(85:105,64:84)],[AS(85:105,85:105)],[AS(85:105,106:126)],[AS(85:105,127:147)],[AS(85:105,148:168)],[AS(85:105,169:189)],[AS(85:105,190:210)];
               [AS(106:126,1:21)],[AS(106:126,22:42)],[AS(106:126,43:63)],[AS(106:126,64:84)],[AS(106:126,85:105)],[AS(106:126,106:126)],[AS(106:126,127:147)],[AS(106:126,148:168)],[AS(106:126,169:189)],[AS(106:126,190:210)];
               [AS(127:147,1:21)],[AS(127:147,22:42)],[AS(127:147,43:63)],[AS(127:147,64:84)],[AS(127:147,85:105)],[AS(127:147,106:126)],[AS(127:147,127:147)],[AS(127:147,148:168)],[AS(127:147,169:189)],[AS(127:147,190:210)];
               [AS(148:168,1:21)],[AS(148:168,22:42)],[AS(148:168,43:63)],[AS(148:168,64:84)],[AS(148:168,85:105)],[AS(148:168,106:126)],[AS(148:168,127:147)],[AS(148:168,148:168)],[AS(148:168,169:189)],[AS(148:168,190:210)];
               [AS(169:189,1:21)],[AS(169:189,22:42)],[AS(169:189,43:63)],[AS(169:189,64:84)],[AS(169:189,85:105)],[AS(169:189,106:126)],[AS(169:189,127:147)],[AS(169:189,148:168)],[AS(169:189,169:189)],[AS(169:189,190:210)];
               [AS(190:210,1:21)],[AS(190:210,22:42)],[AS(190:210,43:63)],[AS(190:210,64:84)],[AS(190:210,85:105)],[AS(190:210,106:126)],[AS(190:210,127:147)],[AS(190:210,148:168)],[AS(190:210,169:189)],[AS(190:210,190:210)]};


 S1 = pop(1:1:n);E1 = pop(n+1:1:2*n);A1 = pop((2*n)+1:1:3*n);J1 = pop((3*n)+1:1:4*n);Q1 = pop((4*n)+1:1:5*n);I1 = pop((5*n)+1:1:6*n);
 H1 = pop((6*n)+1:1:7*n);R1 = pop((7*n)+1:1:8*n);D1 = pop((8*n)+1:1:9*n);Ir1 = pop((9*n)+1:1:10*n);

 S2 = pop((10*n)+1:1:11*n);E2 = pop((11*n)+1:1:12*n);A2 = pop((12*n)+1:1:13*n);J2 = pop((13*n)+1:1:14*n);Q2 = pop((14*n)+1:1:15*n);I2 = pop((15*n)+1:1:16*n);
 H2 = pop((16*n)+1:1:17*n);R2 = pop((17*n)+1:1:18*n);D2 = pop((18*n)+1:1:19*n);Ir2 = pop((19*n)+1:1:20*n);

 S3 = pop((20*n)+1:1:21*n);E3 = pop((21*n)+1:1:22*n);A3 = pop((22*n)+1:1:23*n);J3 = pop((23*n)+1:1:24*n);Q3 = pop((24*n)+1:1:25*n);I3 = pop((25*n)+1:1:26*n);
 H3 = pop((26*n)+1:1:27*n);R3 = pop((27*n)+1:1:28*n);D3 = pop((28*n)+1:1:29*n);Ir3 = pop((29*n)+1:1:30*n);

 S4 = pop((30*n)+1:1:31*n);E4 = pop((31*n)+1:1:32*n);A4 = pop((32*n)+1:1:33*n);J4 = pop((33*n)+1:1:34*n);Q4 = pop((34*n)+1:1:35*n);I4 = pop((35*n)+1:1:36*n);
 H4 = pop((36*n)+1:1:37*n);R4 = pop((37*n)+1:1:38*n);D4 = pop((38*n)+1:1:39*n);Ir4 = pop((39*n)+1:1:40*n);

 S5 = pop((40*n)+1:1:41*n);E5 = pop((41*n)+1:1:42*n);A5 = pop((42*n)+1:1:43*n);J5 = pop((43*n)+1:1:44*n);Q5 = pop((44*n)+1:1:45*n);I5 = pop((45*n)+1:1:46*n);
 H5 = pop((46*n)+1:1:47*n);R5 = pop((47*n)+1:1:48*n);D5 = pop((48*n)+1:1:49*n);Ir5 = pop((49*n)+1:1:50*n);

 S6 = pop((50*n)+1:1:51*n);E6 = pop((51*n)+1:1:52*n);A6 = pop((52*n)+1:1:53*n);J6 = pop((53*n)+1:1:54*n);Q6 = pop((54*n)+1:1:55*n);I6 = pop((55*n)+1:1:56*n);
 H6 = pop((56*n)+1:1:57*n);R6 = pop((57*n)+1:1:58*n);D6 = pop((58*n)+1:1:59*n);Ir6 = pop((59*n)+1:1:60*n);

 S7 = pop((60*n)+1:1:61*n);E7 = pop((61*n)+1:1:62*n);A7 = pop((62*n)+1:1:63*n);J7 = pop((63*n)+1:1:64*n);Q7 = pop((64*n)+1:1:65*n);I7 = pop((65*n)+1:1:66*n);
 H7 = pop((66*n)+1:1:67*n);R7 = pop((67*n)+1:1:68*n);D7 = pop((68*n)+1:1:69*n);Ir7 = pop((69*n)+1:1:70*n);

 S8 = pop((70*n)+1:1:71*n);E8 = pop((71*n)+1:1:72*n);A8 = pop((72*n)+1:1:73*n);J8 = pop((73*n)+1:1:74*n);Q8 = pop((74*n)+1:1:75*n);I8 = pop((75*n)+1:1:76*n);
 H8 = pop((76*n)+1:1:77*n);R8 = pop((77*n)+1:1:78*n);D8 = pop((78*n)+1:1:79*n);Ir8 = pop((79*n)+1:1:80*n);

 S9 = pop((80*n)+1:1:81*n);E9 = pop((81*n)+1:1:82*n);A9 = pop((82*n)+1:1:83*n);J9 = pop((83*n)+1:1:84*n);Q9 = pop((84*n)+1:1:85*n);I9 = pop((85*n)+1:1:86*n);
 H9 = pop((86*n)+1:1:87*n);R9 = pop((87*n)+1:1:88*n);D9 = pop((88*n)+1:1:89*n);Ir9 = pop((89*n)+1:1:90*n);


 S10 = pop((90*n)+1:1:91*n);E10 = pop((91*n)+1:1:92*n);A10 = pop((92*n)+1:1:93*n);J10 = pop((93*n)+1:1:94*n);Q10 = pop((94*n)+1:1:95*n);I10 = pop((95*n)+1:1:96*n);
 H10 = pop((96*n)+1:1:97*n);R10 = pop((97*n)+1:1:98*n);D10 = pop((98*n)+1:1:99*n);Ir10 = pop((99*n)+1:1:100*n);

  %% Store my output here

dpop = zeros(length(pop),1);



%% ODE for each social group

%% SOCIAL GROUP ONE


SocInf1 = ((AgeSoc{1,1}*(J1 + para.iota*Q1 + para.iota*I1 + para.tau*A1)+AgeSoc{1,2}*(J2 + para.iota*Q2 + para.iota*I2 + para.tau*A2)+...
AgeSoc{1,3}*(J3 + para.iota*Q3 + para.iota*I3 + para.tau*A3)+AgeSoc{1,4}*(J4 + para.iota*Q4 + para.iota*I4 + para.tau*A4)+...
AgeSoc{1,5}*(J5 + para.iota*Q5 + para.iota*I5 + para.tau*A5)+AgeSoc{1,6}*(J6 + para.iota*Q6 + para.iota*I6 + para.tau*A6)+...
AgeSoc{1,7}*(J7 + para.iota*Q7 + para.iota*I7 + para.tau*A7)+AgeSoc{1,8}*(J8 + para.iota*Q8 + para.iota*I8 + para.tau*A8)+...
AgeSoc{1,9}*(J9 + para.iota*Q9 + para.iota*I9 + para.tau*A9)+AgeSoc{1,10}*(J10 + para.iota*Q10 + para.iota*I10 + para.tau*A10))).*(S1./para.N1');
 


beta1 =  para.beta*SocInf1;


% Differential Equation

        dpop(1:1:n) = - (beta1)+ para.epsilon*R1;
        dpop(n+1:1:2*n) =   (beta1)- pi(m).*E1 - (1-pi(m))*para.theta*E1 - (1-pi(m))*para.rho*E1;
        dpop((2*n)+1:1:3*n) = (1-pi(m))*para.rho*E1 - para.gamma*A1;
        dpop((3*n)+1:1:4*n) = (1-pi(m))*para.theta*E1 - para.delta.*J1 - para.gamma*J1;
        dpop((4*n)+1:1:5*n) = pi(m)*E1 - 1/para.PHI*para.p*Q1 - 1/para.PHI*(1-para.p)*para.gamma*Q1;
        dpop((5*n)+1:1:6*n) = 1/para.PHI*para.p*Q1 - para.delta.*I1 -para.gamma*I1;
        dpop((6*n)+1:1:7*n) = para.delta.*J1+ para.delta.*I1 - para.d.*H1 - para.gamma_2*H1;
        dpop((7*n)+1:1:8*n) = para.gamma*J1 + para.gamma*A1 +para.gamma_2*H1 +1/para.PHI*(1-para.p)*para.gamma*Q1 + para.gamma*I1 - para.epsilon*R1;
        dpop((8*n)+1:1:9*n)=  para.d.*H1;
        dpop((9*n)+1:1:10*n) = 1/para.PHI*para.p*Q1;

        %% SOCIAL GROUP TWO
        

SocInf2 = (AgeSoc{2,1}*(J1 + para.iota*Q1 + para.iota*I1 + para.tau*A1)+AgeSoc{2,2}*(J2 + para.iota*Q2 + para.iota*I2 + para.tau*A2)+...
AgeSoc{2,3}*(J3 + para.iota*Q3 + para.iota*I3 + para.tau*A3)+AgeSoc{2,4}*(J4 + para.iota*Q4 + para.iota*I4 + para.tau*A4)+...
AgeSoc{2,5}*(J5 + para.iota*Q5 + para.iota*I5 + para.tau*A5)+AgeSoc{2,6}*(J6 + para.iota*Q6 + para.iota*I6 + para.tau*A6)+...
AgeSoc{2,7}*(J7 + para.iota*Q7 + para.iota*I7 + para.tau*A7)+AgeSoc{2,8}*(J8 + para.iota*Q8 + para.iota*I8 + para.tau*A8)+...
AgeSoc{2,9}*(J9 + para.iota*Q9 + para.iota*I9 + para.tau*A9)+AgeSoc{2,10}*(J10 + para.iota*Q10 + para.iota*I10 + para.tau*A10)).* (S2./para.N2');

% AgeMat'* sum(w(2,:));


beta2 =  para.beta*SocInf2;


% Differential Equation

        dpop((10*n)+1:1:11*n) = - beta2 + para.epsilon*R2;
        dpop((11*n)+1:1:12*n) =  beta2 - pi(m).*E2 - (1-pi(m))*para.theta*E2 - (1-pi(m))*para.rho*E2;
        dpop((12*n)+1:1:13*n) = (1-pi(m))*para.rho*E2 - para.gamma*A2;
        dpop((13*n)+1:1:14*n) = (1-pi(m))*para.theta*E2 - para.delta.*J2 - para.gamma*J2;
        dpop((14*n)+1:1:15*n) = pi(m)*E2 - 1/para.PHI*para.p*Q2 - 1/para.PHI*(1-para.p)*para.gamma*Q2;
        dpop((15*n)+1:1:16*n) = 1/para.PHI*para.p*Q2 - para.delta.*I2 -para.gamma*I2;
        dpop((16*n)+1:1:17*n) = para.delta.*J2+ para.delta.*I2 - para.d.*H2 - para.gamma_2*H2;
        dpop((17*n)+1:1:18*n) = para.gamma*J2 + para.gamma*A2 +para.gamma_2*H2 +1/para.PHI*(1-para.p)*para.gamma*Q2 + para.gamma*I2 - para.epsilon*R2;
        dpop((18*n)+1:1:19*n)=  para.d.*H2;
        dpop((19*n)+1:1:20*n) = 1/para.PHI*para.p*Q2;

 %% SOCIAL GROUP THREE

SocInf3 = (AgeSoc{3,1}*(J1 + para.iota*Q1 + para.iota*I1 + para.tau*A1)+AgeSoc{3,2}*(J2 + para.iota*Q2 + para.iota*I2 + para.tau*A2)+...
AgeSoc{3,3}*(J3 + para.iota*Q3 + para.iota*I3 + para.tau*A3)+AgeSoc{3,4}*(J4 + para.iota*Q4 + para.iota*I4 + para.tau*A4)+...
AgeSoc{3,5}*(J5 + para.iota*Q5 + para.iota*I5 + para.tau*A5)+AgeSoc{3,6}*(J6 + para.iota*Q6 + para.iota*I6 + para.tau*A6)+...
AgeSoc{3,7}*(J7 + para.iota*Q7 + para.iota*I7 + para.tau*A7)+AgeSoc{3,8}*(J8 + para.iota*Q8 + para.iota*I8 + para.tau*A8)+...
AgeSoc{3,9}*(J9 + para.iota*Q9 + para.iota*I9 + para.tau*A9)+AgeSoc{3,10}*(J10 + para.iota*Q10 + para.iota*I10 + para.tau*A10)).* (S3./para.N3');

% AgeMat'* sum(w(3,:));



beta3 =   para.beta*SocInf3;
% SocInf3*(S3./para.N3').*(para.tau*J3 + para.iota*Q3 + para.iota*I3 + A3);

% Differential Equation

        dpop((20*n)+1:1:21*n) =  -beta3+para.epsilon * R3;
        dpop((21*n)+1:1:22*n) =   beta3 - pi(m).*E3 - (1-pi(m))*para.theta*E3 - (1-pi(m))*para.rho*E3;
        dpop((22*n)+1:1:23*n) = (1-pi(m))*para.rho*E3 - para.gamma*A3;
        dpop((23*n)+1:1:24*n) = (1-pi(m))*para.theta*E3 - para.delta.*J3 - para.gamma*J3;
        dpop((24*n)+1:1:25*n) = pi(m)*E3 - 1/para.PHI*para.p*Q3 - 1/para.PHI*(1-para.p)*para.gamma*Q3;
        dpop((25*n)+1:1:26*n) = 1/para.PHI*para.p*Q3 - para.delta.*I3 -para.gamma*I3;
        dpop((26*n)+1:1:27*n) = para.delta.*J3+ para.delta.*I3 - para.d.*H3 - para.gamma_2*H3;
        dpop((27*n)+1:1:28*n) = para.gamma*J3 + para.gamma*A3 +para.gamma_2*H3 +1/para.PHI*(1-para.p)*para.gamma*Q3 + para.gamma*I3 - para.epsilon*R3;
        dpop((28*n)+1:1:29*n)=  para.d.*H3;
        dpop((29*n)+1:1:30*n) = 1/para.PHI*para.p*Q3;

         %% SOCIAL GROUP FOUR

SocInf4 = (AgeSoc{4,1}*(J1 + para.iota*Q1 + para.iota*I1 + para.tau*A1)+AgeSoc{4,2}*(J2 + para.iota*Q2 + para.iota*I2 + para.tau*A2)+...
AgeSoc{4,3}*(J3 + para.iota*Q3 + para.iota*I3 + para.tau*A3)+AgeSoc{4,4}*(J4 + para.iota*Q4 + para.iota*I4 + para.tau*A4)+...
AgeSoc{4,5}*(J5 + para.iota*Q5 + para.iota*I5 + para.tau*A5)+AgeSoc{4,6}*(J6 + para.iota*Q6 + para.iota*I6 + para.tau*A6)+...
AgeSoc{4,7}*(J7 + para.iota*Q7 + para.iota*I7 + para.tau*A7)+AgeSoc{4,8}*(J8 + para.iota*Q8 + para.iota*I8 + para.tau*A8)+...
AgeSoc{4,9}*(J9 + para.iota*Q9 + para.iota*I9 + para.tau*A9)+AgeSoc{4,10}*(J10 + para.iota*Q10 + para.iota*I10 + para.tau*A10)).* (S4./para.N4');


beta4 =  para.beta*SocInf4;


% Differential Equation

        dpop((30*n)+1:1:31*n) =   -beta4+para.epsilon * R4;
        dpop((31*n)+1:1:32*n) =   beta4 - pi(m).*E4 - (1-pi(m))*para.theta*E4 - (1-pi(m))*para.rho*E4;
        dpop((32*n)+1:1:33*n) = (1-pi(m))*para.rho*E4 - para.gamma*A4;
        dpop((33*n)+1:1:34*n) = (1-pi(m))*para.theta*E4 - para.delta.*J4 - para.gamma*J4;
        dpop((34*n)+1:1:35*n) = pi(m)*E4 - 1/para.PHI*para.p*Q4 - 1/para.PHI*(1-para.p)*para.gamma*Q4;
        dpop((35*n)+1:1:36*n) = 1/para.PHI*para.p*Q4 - para.delta.*I4 -para.gamma*I4;
        dpop((36*n)+1:1:37*n) = para.delta.*J4+ para.delta.*I4 - para.d.*H4 - para.gamma_2*H4;
        dpop((37*n)+1:1:38*n) = para.gamma*J4 + para.gamma*A4 +para.gamma_2*H4 +1/para.PHI*(1-para.p)*para.gamma*Q4 + para.gamma*I4 - para.epsilon*R4;
        dpop((38*n)+1:1:39*n) =  para.d.*H4;
        dpop((39*n)+1:1:40*n) = 1/para.PHI*para.p*Q4;


         %% SOCIAL GROUP FIVE

SocInf5 = (AgeSoc{5,1}*(J1 + para.iota*Q1 + para.iota*I1 + para.tau*A1)+AgeSoc{5,2}*(J2 + para.iota*Q2 + para.iota*I2 + para.tau*A2)+...
AgeSoc{5,3}*(J3 + para.iota*Q3 + para.iota*I3 + para.tau*A3)+AgeSoc{5,4}*(J4 + para.iota*Q4 + para.iota*I4 + para.tau*A4)+...
AgeSoc{5,5}*(J5 + para.iota*Q5 + para.iota*I5 + para.tau*A5)+AgeSoc{5,6}*(J6 + para.iota*Q6 + para.iota*I6 + para.tau*A6)+...
AgeSoc{5,7}*(J7 + para.iota*Q7 + para.iota*I7 + para.tau*A7)+AgeSoc{5,8}*(J8 + para.iota*Q8 + para.iota*I8 + para.tau*A8)+...
AgeSoc{5,9}*(J9 + para.iota*Q9 + para.iota*I9 + para.tau*A9)+AgeSoc{5,10}*(J10 + para.iota*Q10 + para.iota*I10 + para.tau*A10)).* (S5./para.N5');



beta5 = para.beta*SocInf5;

% Differential Equation

         dpop((40*n)+1:1:41*n) =  -beta5 +para.epsilon * R5;
        dpop((41*n)+1:1:42*n) =   beta5  - pi(m).*E5 - (1-pi(m))*para.theta*E5 - (1-pi(m))*para.rho*E5;
        dpop((42*n)+1:1:43*n) = (1-pi(m))*para.rho*E5 - para.gamma*A5;
        dpop((43*n)+1:1:44*n) = (1-pi(m))*para.theta*E5 - para.delta.*J5 - para.gamma*J5;
        dpop((44*n)+1:1:45*n) = pi(m)*E5 - 1/para.PHI*para.p*Q5 - 1/para.PHI*(1-para.p)*para.gamma*Q5;
        dpop((45*n)+1:1:46*n) = 1/para.PHI*para.p*Q5 - para.delta.*I5 -para.gamma*I5;
        dpop((46*n)+1:1:47*n) = para.delta.*J5+ para.delta.*I5 - para.d.*H5 - para.gamma_2*H5;
        dpop((47*n)+1:1:48*n) = para.gamma*J5 + para.gamma*A5 +para.gamma_2*H5 +1/para.PHI*(1-para.p)*para.gamma*Q5 + para.gamma*I5 - para.epsilon*R5;
        dpop((48*n)+1:1:49*n)=  para.d.*H5;
        dpop((49*n)+1:1:50*n) = 1/para.PHI*para.p*Q5;


         %% SOCIAL GROUP SIX

SocInf6 = (AgeSoc{6,1}*(J1 + para.iota*Q1 + para.iota*I1 + para.tau*A1)+AgeSoc{6,2}*(J2 + para.iota*Q2 + para.iota*I2 + para.tau*A2)+...
AgeSoc{6,3}*(J3 + para.iota*Q3 + para.iota*I3 + para.tau*A3)+AgeSoc{6,4}*(J4 + para.iota*Q4 + para.iota*I4 + para.tau*A4)+...
AgeSoc{6,5}*(J5 + para.iota*Q5 + para.iota*I5 + para.tau*A5)+AgeSoc{6,6}*(J6 + para.iota*Q6 + para.iota*I6 + para.tau*A6)+...
AgeSoc{6,7}*(J7 + para.iota*Q7 + para.iota*I7 + para.tau*A7)+AgeSoc{6,8}*(J8 + para.iota*Q8 + para.iota*I8 + para.tau*A8)+...
AgeSoc{6,9}*(J9 + para.iota*Q9 + para.iota*I9 + para.tau*A9)+AgeSoc{6,10}*(J10 + para.iota*Q10 + para.iota*I10 + para.tau*A10)).* (S6./para.N6');


beta6 =  para.beta*SocInf6;


% Differential Equation

        dpop((50*n)+1:1:51*n) =  -beta6 +para.epsilon * R6;
        dpop((51*n)+1:1:52*n) =   beta6  - pi(m).*E6 - (1-pi(m))*para.theta*E6 - (1-pi(m))*para.rho*E6;
        dpop((52*n)+1:1:53*n) = (1-pi(m))*para.rho*E6 - para.gamma*A6;
        dpop((53*n)+1:1:54*n) = (1-pi(m))*para.theta*E6 - para.delta.*J6 - para.gamma*J6;
        dpop((54*n)+1:1:55*n) = pi(m)*E6 - 1/para.PHI*para.p*Q6 - 1/para.PHI*(1-para.p)*para.gamma*Q6;
        dpop((55*n)+1:1:56*n) = 1/para.PHI*para.p*Q6 - para.delta.*I6 -para.gamma*I6;
        dpop((56*n)+1:1:57*n) = para.delta.*J6+ para.delta.*I6 - para.d.*H6 - para.gamma_2*H6;
        dpop((57*n)+1:1:58*n) = para.gamma*J6 + para.gamma*A6 +para.gamma_2*H6 +1/para.PHI*(1-para.p)*para.gamma*Q6 + para.gamma*I6 - para.epsilon*R6;
        dpop((58*n)+1:1:59*n)=  para.d.*H6;
        dpop((59*n)+1:1:60*n) = 1/para.PHI*para.p*Q6;


         %% SOCIAL GROUP SEVEN

SocInf7 = (AgeSoc{7,1}*(J1 + para.iota*Q1 + para.iota*I1 + para.tau*A1)+AgeSoc{7,2}*(J2 + para.iota*Q2 + para.iota*I2 + para.tau*A2)+...
AgeSoc{7,3}*(J3 + para.iota*Q3 + para.iota*I3 + para.tau*A3)+AgeSoc{7,4}*(J4 + para.iota*Q4 + para.iota*I4 + para.tau*A4)+...
AgeSoc{7,5}*(J5 + para.iota*Q5 + para.iota*I5 + para.tau*A5)+AgeSoc{7,6}*(J6 + para.iota*Q6 + para.iota*I6 + para.tau*A6)+...
AgeSoc{7,7}*(J7 + para.iota*Q7 + para.iota*I7 + para.tau*A7)+AgeSoc{7,8}*(J8 + para.iota*Q8 + para.iota*I8 + para.tau*A8)+...
AgeSoc{7,9}*(J9 + para.iota*Q9 + para.iota*I9 + para.tau*A9)+AgeSoc{7,10}*(J10 + para.iota*Q10 + para.iota*I10 + para.tau*A10)).* (S7./para.N7');


% 
beta7 =   para.beta*SocInf7;

% Differential Equation

        dpop((60*n)+1:1:61*n) =  -beta7 +para.epsilon * R7;
        dpop((61*n)+1:1:62*n) =   beta7 - pi(m).*E7 - (1-pi(m))*para.theta*E7 - (1-pi(m))*para.rho*E7;
        dpop((62*n)+1:1:63*n) = (1-pi(m))*para.rho*E7 - para.gamma*A7;
        dpop((63*n)+1:1:64*n) = (1-pi(m))*para.theta*E7 - para.delta.*J7 - para.gamma*J7;
        dpop((64*n)+1:1:65*n) = pi(m)*E7 - 1/para.PHI*para.p*Q7 - 1/para.PHI*(1-para.p)*para.gamma*Q7;
        dpop((65*n)+1:1:66*n) = 1/para.PHI*para.p*Q7 - para.delta.*I7 -para.gamma*I7;
        dpop((66*n)+1:1:67*n) = para.delta.*J7+ para.delta.*I7 - para.d.*H7 - para.gamma_2*H7;
        dpop((67*n)+1:1:68*n) = para.gamma*J7 + para.gamma*A7 +para.gamma_2*H7 +1/para.PHI*(1-para.p)*para.gamma*Q7 + para.gamma*I7 - para.epsilon*R7;
        dpop((68*n)+1:1:69*n)=  para.d.*H7;
        dpop((69*n)+1:1:70*n) = 1/para.PHI*para.p*Q7;


         %% SOCIAL GROUP EIGHT

SocInf8 = (AgeSoc{8,1}*(J1 + para.iota*Q1 + para.iota*I1 + para.tau*A1)+AgeSoc{8,2}*(J2 + para.iota*Q2 + para.iota*I2 + para.tau*A2)+...
AgeSoc{8,3}*(J3 + para.iota*Q3 + para.iota*I3 + para.tau*A3)+AgeSoc{8,4}*(J4 + para.iota*Q4 + para.iota*I4 + para.tau*A4)+...
AgeSoc{8,5}*(J5 + para.iota*Q5 + para.iota*I5 + para.tau*A5)+AgeSoc{8,6}*(J6 + para.iota*Q6 + para.iota*I6 + para.tau*A6)+...
AgeSoc{8,7}*(J7 + para.iota*Q7 + para.iota*I7 + para.tau*A7)+AgeSoc{8,8}*(J8 + para.iota*Q8 + para.iota*I8 + para.tau*A8)+...
AgeSoc{8,9}*(J9 + para.iota*Q9 + para.iota*I9 + para.tau*A9)+AgeSoc{8,10}*(J10 + para.iota*Q10 + para.iota*I10 + para.tau*A10)).* (S8./para.N8');


beta8 =   para.beta*SocInf8;

% Differential Equation
 
        dpop((70*n)+1:1:71*n) =   -beta8 +para.epsilon * R8;
        dpop((71*n)+1:1:72*n) =   beta8  - pi(m).*E8 - (1-pi(m))*para.theta*E8 - (1-pi(m))*para.rho*E8;
        dpop((72*n)+1:1:73*n) = (1-pi(m))*para.rho*E8 - para.gamma*A8;
        dpop((73*n)+1:1:74*n) = (1-pi(m))*para.theta*E8 - para.delta.*J8 - para.gamma*J8;
        dpop((74*n)+1:1:75*n) = pi(m)*E8 - 1/para.PHI*para.p*Q8 - 1/para.PHI*(1-para.p)*para.gamma*Q8;
        dpop((75*n)+1:1:76*n) = 1/para.PHI*para.p*Q8 - para.delta.*I8 -para.gamma*I8;
        dpop((76*n)+1:1:77*n) = para.delta.*J8+ para.delta.*I8 - para.d.*H8 - para.gamma_2*H8;
        dpop((77*n)+1:1:78*n) = para.gamma*J8 + para.gamma*A8 +para.gamma_2*H8 +1/para.PHI*(1-para.p)*para.gamma*Q8 + para.gamma*I8 - para.epsilon*R8;
        dpop((78*n)+1:1:79*n) =  para.d.*H8;
        dpop((79*n)+1:1:80*n) = 1/para.PHI*para.p*Q8;


         %% SOCIAL GROUP NINE

SocInf9 = (AgeSoc{9,1}*(J1 + para.iota*Q1 + para.iota*I1 + para.tau*A1)+AgeSoc{9,2}*(J2 + para.iota*Q2 + para.iota*I2 + para.tau*A2)+...
AgeSoc{9,3}*(J3 + para.iota*Q3 + para.iota*I3 + para.tau*A3)+AgeSoc{9,4}*(J4 + para.iota*Q4 + para.iota*I4 + para.tau*A4)+...
AgeSoc{9,5}*(J5 + para.iota*Q5 + para.iota*I5 + para.tau*A5)+AgeSoc{9,6}*(J6 + para.iota*Q6 + para.iota*I6 + para.tau*A6)+...
AgeSoc{9,7}*(J7 + para.iota*Q7 + para.iota*I7 + para.tau*A7)+AgeSoc{9,8}*(J8 + para.iota*Q8 + para.iota*I8 + para.tau*A8)+...
AgeSoc{9,9}*(J9 + para.iota*Q9 + para.iota*I9 + para.tau*A9)+AgeSoc{9,10}*(J10 + para.iota*Q10 + para.iota*I10 + para.tau*A10)).* (S9./para.N9');


beta9 =  para.beta*SocInf9;

% Differential Equation


        dpop((80*n)+1:1:81*n) =  -beta9 +para.epsilon * R9;
        dpop((81*n)+1:1:82*n) =   beta9 - pi(m).*E9 - (1-pi(m))*para.theta*E9 - (1-pi(m))*para.rho*E9;
        dpop((82*n)+1:1:83*n) = (1-pi(m))*para.rho*E9 - para.gamma*A9;
        dpop((83*n)+1:1:84*n) = (1-pi(m))*para.theta*E9 - para.delta.*J9 - para.gamma*J9;
        dpop((84*n)+1:1:85*n) = pi(m)*E9 - 1/para.PHI*para.p*Q9 - 1/para.PHI*(1-para.p)*para.gamma*Q9;
        dpop((85*n)+1:1:86*n) = 1/para.PHI*para.p*Q9 - para.delta.*I9 -para.gamma*I9;
        dpop((86*n)+1:1:87*n) = para.delta.*J9+ para.delta.*I9 - para.d.*H9 - para.gamma_2*H9;
        dpop((87*n)+1:1:88*n) = para.gamma*J9 + para.gamma*A9 +para.gamma_2*H9 +1/para.PHI*(1-para.p)*para.gamma*Q9 + para.gamma*I9 - para.epsilon*R9;
        dpop((88*n)+1:1:89*n)=  para.d.*H9;
        dpop((89*n)+1:1:90*n) = 1/para.PHI*para.p*Q9;


         %% SOCIAL GROUP TEN

SocInf10 = (AgeSoc{10,1}*(J1 + para.iota*Q1 + para.iota*I1 + para.tau*A1)+AgeSoc{10,2}*(J2 + para.iota*Q2 + para.iota*I2 + para.tau*A2)+...
            AgeSoc{10,3}*(J3 + para.iota*Q3 + para.iota*I3 + para.tau*A3)+AgeSoc{10,4}*(J4 + para.iota*Q4 + para.iota*I4 + para.tau*A4)+...
            AgeSoc{10,5}*(J5 + para.iota*Q5 + para.iota*I5 + para.tau*A5)+AgeSoc{10,6}*(J6 + para.iota*Q6 + para.iota*I6 + para.tau*A6)+...
            AgeSoc{10,7}*(J7 + para.iota*Q7 + para.iota*I7 + para.tau*A7)+AgeSoc{10,8}*(J8 + para.iota*Q8 + para.iota*I8 + para.tau*A8)+...
            AgeSoc{10,9}*(J9 + para.iota*Q9 + para.iota*I9 + para.tau*A9)+AgeSoc{10,10}*(J10 + para.iota*Q10 + para.iota*I10 + para.tau*A10)).* (S10./para.N10');


beta10 =  para.beta*SocInf10;

% Differential Equation


        dpop((90*n)+1:1:91*n) =  -beta10 +para.epsilon * R10;
        dpop((91*n)+1:1:92*n) =   beta10 - pi(m).*E10 - (1-pi(m))*para.theta*E10 - (1-pi(m))*para.rho*E10;
        dpop((92*n)+1:1:93*n) = (1-pi(m))*para.rho*E10 - para.gamma*A10;
        dpop((93*n)+1:1:94*n) = (1-pi(m))*para.theta*E10 - para.delta.*J10 - para.gamma*J10;
        dpop((94*n)+1:1:95*n) = pi(m)*E10 - 1/para.PHI*para.p*Q10 - 1/para.PHI*(1-para.p)*para.gamma*Q10;
        dpop((95*n)+1:1:96*n) = 1/para.PHI*para.p*Q10 - para.delta.*I10 -para.gamma*I10;
        dpop((96*n)+1:1:97*n) = para.delta.*J10+ para.delta.*I10 - para.d.*H10 - para.gamma_2*H10;
        dpop((97*n)+1:1:98*n) = para.gamma*J10 + para.gamma*A10 +para.gamma_2*H10 +1/para.PHI*(1-para.p)*para.gamma*Q10 + para.gamma*I10 - para.epsilon*R10;
        dpop((98*n)+1:1:99*n) =  para.d.*H10;
        dpop((99*n)+1:1:100*n) = 1/para.PHI*para.p*Q10;
 end
end


