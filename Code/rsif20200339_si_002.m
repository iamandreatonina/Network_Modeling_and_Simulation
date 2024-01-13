% Mathematical modeling of breast cancer cells in response to endocrine therapy and Cdk4/6 inhibition
% Journal of the Royal Society Interface
% Wei He, Diane M. Demas, Isabel P. Conde, Ayesha N. Shajahan-Haq, William T. Baumann

% The main code runs the simulation. It will give the simulation results for all
% proteins at -E2, E2+ICI(500nM), E2+palbo(1uM), -E2+ICI(500nM) and
% -E2+palbo(1uM) conditions and the cell number at E2 control, -E2,
% E2+palbo(1uM), -E2+ICI(500nM) and -E2+palbo(1uM)
clc;clear;
% Protein values from western blot
OutE2depprotein.mean = [1,0.750,0.875,0.828,0.209,0.188;...
    1,0.982,0.945,1.044,0.860,0.921;...
    1,1.141,2.03,2.371,2.342,2.692;...
    1,1.0156,1.134,1.163,1.249,1.337;...
    1,0.976,0.892,0.844,0.573,0.547;...
    1,1.218,0.433,0.729,0.387,0.413];
OutE2depprotein.sem = [0,0.108,0.130,0.057,0.007,0.018;...
    0,0.110,0.131,0.169,0.124,0.072;...
    0,0.174,0.324,0.465,0.382,0.238;...
    0,0.286,0.422,0.32,0.337,0.344;...
    0,0.016,0.055,0.079,0.057,0.119;...
    0,0.154,0.049,0.091,0.083,0.200];
OutE2ICI500nMprotein.mean = [1,0.645,0.482,0.350,0.183,0.158;...
    1,0.946,0.820,0.845,0.761,0.775;...
    1,0.693,0.675,0.629,0.691,0.720;...
    1,1.862,1.837,1.738,1.580,1.520;...
    1,0.952,0.836,0.663,0.603,0.727;...
    1,1.167,0.373,0.099,0.314,0.294];
OutE2ICI500nMprotein.sem = [0,0.039,0.076,0.054,0.021,0.054;...
    0,0.123,0.155,0.124,0.125,0.100;...
    0,0.087,0.052,0.024,0.004,0.045;...
    0,0.628,0.473,0.893,0.846,0.811;...
    0,0.027,0.115,0.123,0.095,0.143;...
    0,0.304,0.119,0.054,0.141,0.117];
OutE2depICI500nMprotein.mean = [1,0.562,0.355,0.242,0.118,0.104;...
    1,1.082,0.862,0.789,0.653,0.678;...
    1,0.681,0.593,0.510,0.470,0.522;...
    1,1.251,1.618,1.292,1.042,1.256;...
    1,0.937,0.528,0.370,0.322,0.317;...
    1,0.963,0.300,0.119,0.048,0.057];
OutE2depICI500nMprotein.sem = [0,0.025,0.056,0.019,0.010,0.017;...
    0,0.109,0.161,0.146,0.072,0.115;...
    0,0.016,0.035,0.070,0.053,0.076;...
    0,0.208,0.172,0.286,0.255,0.455;...
    0,0.031,0.107,0.080,0.056,0.100;...
    0,0.066,0.107,0.062,0.018,0.038];
OutE2palbo1uMprotein.mean = [1,0.922,0.585,0.379,0.330;...
    1,0.953,0.808,0.553,0.451];
OutE2palbo1uMprotein.sem = [0,0.071,0.013,0.024,0.018;...
    0,0.018,0.061,0.052,0.020];
OutE2deppalbo1uMprotein.mean = [1,0.613,0.346,0.101,0.089;...
    1,0.634,0.489,0.107,0.140];
OutE2deppalbo1uMprotein.sem = [0,0.047,0.043,0.013,0.013;...
    0,0.058,0.023,0.0225,0.020];
% Cell numbers from coulter counter
OutE2proCounter.mean = [1,1.748,5.372,19.347,22.456];
OutE2proCounter.sem = [0,0.048,0.306,1.630,1.318];
OutE2depproCounter.mean = [1,1.322,3.772,8.965,9.938];
OutE2depproCounter.sem = [0,0.049,0.670,1.309,1.199];
OutE2ICI500nMproCounter.mean = [1,1.344,2.524,3.579,3.840];
OutE2ICI500nMproCounter.sem = [0,0.135,0.169,0.241,0.241];
OutE2depICI500nMproCounter.mean = [1,1.269,2.325,2.945,3.103];
OutE2depICI500nMproCounter.sem = [0,0.078,0.267,0.165,0.117];
OutE2palbo1uMproCounter.mean = [1,1.569,2.558,3.362,3.421];
OutE2palbo1uMproCounter.sem = [0,0.052,0.264,0.576,0.414];
OutE2deppalbo1uMproCounter.mean = [1,1.144,1.460,1.230,1.427];
OutE2deppalbo1uMproCounter.sem = [0,0.081,0.135,0.191,0.136];

% protein name
proteinname{1} = 'c-Myc';proteinname{2} = 'total cyclinD1';
proteinname{3} = 'ER'; proteinname{4} = 'p21';
proteinname{5} = 'RB1-pp';proteinname{6} = 'total RB1';

proteinnamepalbo{1} = 'c-Myc';proteinnamepalbo{2} = 'RB1-pp';
ttimedata = [0, 0.1667, 1, 3, 5, 7]; % day
tpro = [0,1,3,6,7];
Numvar = 14;
ValE2normal = 0.01; % E2 concentration 0.01uM
ValICI100nM = 0.5; % ICI concentration 0.5uM
Valpalbo1000nM = 1; % Palbo concentration 1uM
ind_E2dep = 58;
odefun = @ode23tb;
E2dep = 0;

% initial value of variables
x0 = ones(Numvar,1);
x0(3) = 0; % ICIER = 0
x0(6) = 0; % cyclinD1palbo = 0 
x0(end) = 1;
% Parameter used
PAR = [0.0207;0.1;0.3;4266.4776;1;206.8202;1;0.5184;0.1375;1.73;11.5734;2.6086;0.1704;32.0864;1;5.1776;...
    1;0.33;2.31;14.1779;0.0282;6.998;37.3017;5.5045;5.3598;1.2879;1.39;9.2431;2.4352;0.8058;1.39;20.9396;...
    1;3.267;0.35;390.8695;10.7606;7;15.0134;23.6192;0.8;5.3707;9.2727;0.05;14.0135;3.1884;0.6772;7;6.3544;...
    0.3503;59.8461;2.11;0.1196;266.3026;6.0056;6;37.3859;9.6676e-05;2.3e-05];

n = 1:size(OutE2ICI500nMprotein.mean,1);

E2 = ValE2normal; % concentration of E2
ICI = 0; % value of ICI, if bigger than 0, means add ICI
palbo = 0; % concentration of palbo

tspan = [0,1e3];
opt = odeset('RelTol',1e-6,'AbsTol',1e-6);
[t,x] = odefun(@(t,X)E2ICIpalbo_model(t,X,PAR,ICI,E2,palbo,E2dep),tspan,x0,opt);
% set the steady state
x0 = x(end,:);
%% E2 control
clc;close
tspan = 0:200;
x0(3) = 0;x0(6) = 0;x0(end) = 1;
[tE2,xE2] = odefun(@(t,X)E2ICIpalbo_model(t,X,PAR,ICI,E2,palbo,E2dep),tspan,x0,opt);
simresultE2 = E2ICIpalbo_assignval(xE2);

% +E2+ICI(500nM)
E2 = ValE2normal; % concentration of E2
ICI = ValICI100nM; % concentration of ICI
palbo = 0; % concentration of palbo
[tE2ICI,xE2ICI] = odefun(@(t,X)E2ICIpalbo_model(t,X,PAR,ICI,E2,palbo,E2dep),tspan,x0,opt);
simresultE2ICI = E2ICIpalbo_assignval(xE2ICI);

% -E2
E2 = PAR(ind_E2dep); % concentration of E2 in deprived medium
ICI = 0; % ICI concentration
palbo = 0; % concentration of palbo
E2dep = 1;
[tE2dep,xE2dep] = odefun(@(t,X)E2ICIpalbo_model(t,X,PAR,ICI,E2,palbo,E2dep),tspan,x0,opt);
simresultE2dep = E2ICIpalbo_assignval(xE2dep);

% -E2+ICI(500nM)
E2 = PAR(ind_E2dep); % E2 deprived
ICI = ValICI100nM; % ICI concentration
palbo = 0; % concentration of palbo
[tE2depICI,xE2depICI] = odefun(@(t,X)E2ICIpalbo_model(t,X,PAR,ICI,E2,palbo,E2dep),tspan,x0,opt);
simresultE2depICI = E2ICIpalbo_assignval(xE2depICI);

% +E2+palbo(1uM)
E2 = ValE2normal; % concentration of E2
ICI = 0; % concentration of ICI
palbo = Valpalbo1000nM; % concentration of palbo
E2dep = 0;
[tE2palbo,xE2palbo] = odefun(@(t,X)E2ICIpalbo_model(t,X,PAR,ICI,E2,palbo,E2dep),tspan,x0,opt);
simresultE2palbo = E2ICIpalbo_assignval(xE2palbo);

% -E2+palbo(1uM)
E2 = PAR(ind_E2dep); % concentration of E2 in deprived medium
ICI = 0; % concentration of ICI
E2dep = 1;
palbo = Valpalbo1000nM; % concentration of palbo
[tE2deppalbo,xE2deppalbo] = odefun(@(t,X)E2ICIpalbo_model(t,X,PAR,ICI,E2,palbo,E2dep),tspan,x0,opt);
simresultE2deppalbo = E2ICIpalbo_assignval(xE2deppalbo);
%% protein data
% plots of all proteins of -E2
for i = n
    figure;
    errorbar(ttimedata,OutE2depprotein.mean(i,:),OutE2depprotein.sem(i,:),'-ro','linewidth',2);hold on
    plot(tE2dep/24,simresultE2dep(i,:),'c','linewidth',2)
    title(['-E2 ',proteinname{i}])
    set(gca,'fontsize',20,'linewidth',2,'xtick',[0 0.1667 1.0000 3.0000 5.0000 7.0000],...
        'xticklabel',{'','4h','1d','3d','5d','7d'})
    legend({'Experiment','Simulation'},'location','NW','box','off')
    xlim([-0.2,7.2])
    ylim([0,4])
end
% plots of all proteins of +E2+ICI(500nM)
for i = n
    figure;
    errorbar(ttimedata,OutE2ICI500nMprotein.mean(i,:),OutE2ICI500nMprotein.sem(i,:),'-ro','linewidth',2);hold on
    plot(tE2ICI/24,simresultE2ICI(i,:),'c','linewidth',2)
    title(['+E2+ICI(500nM) ',proteinname{i}])
    set(gca,'fontsize',20,'linewidth',2,'xtick',[0 0.1667 1.0000 3.0000 5.0000 7.0000],...
        'xticklabel',{'','4h','1d','3d','5d','7d'})
    legend({'Experiment','Simulation'},'location','NW','box','off')
    xlim([-0.2,7.2])
    ylim([0,4])
end

% plots of all proteins of -E2+ICI(500nM)
for i = n
    figure
    errorbar(ttimedata,OutE2depICI500nMprotein.mean(i,:),OutE2depICI500nMprotein.sem(i,:),'-ro','linewidth',2);hold on
    plot(tE2depICI/24,simresultE2depICI(i,:),'b','linewidth',2)
    title(['-E2+ICI(500nM) ',proteinname{i},])
    set(gca,'fontsize',20,'linewidth',2,'xtick',[0 0.1667 1.0000 3.0000 5.0000 7.0000],...
        'xticklabel',{'','4h','1d','3d','5d','7d'})
    legend({'Experiment','Prediction'},'location','NW','box','off')
    xlim([-0.2,7.2])
    ylim([0,2.5])
end
% plots of all proteins of +E2+palbo(1uM)
simresultE2palboplot = simresultE2palbo([1,5],:);
for i = 1 : size(simresultE2palboplot,1)
    figure
    errorbar(tpro,OutE2palbo1uMprotein.mean(i,:),OutE2palbo1uMprotein.sem(i,:),'-ro','linewidth',2);hold on
    plot(tE2palbo/24,simresultE2palboplot(i,:),'c','linewidth',2)
    title(['+E2+palbo(1uM) ',proteinnamepalbo{i},])
    set(gca,'fontsize',20,'linewidth',2,'xtick',[0 1.0000 3.0000 6.0000 7.0000],...
        'xticklabel',{'0h','1d','3d','6d','7d'})
    legend({'Experiment','Simulation'},'location','NW','box','off')
    xlim([-0.2,7.2])
    ylim([0,2.5])
end
% plots of all proteins of -E2+palbo(1uM)
simresultE2deppalboplot = simresultE2deppalbo([1,5],:);
for i = 1 : size(simresultE2deppalboplot,1)
    figure
    errorbar(tpro,OutE2deppalbo1uMprotein.mean(i,:),OutE2deppalbo1uMprotein.sem(i,:),'-ro','linewidth',2);hold on
    plot(tE2deppalbo/24,simresultE2deppalboplot(i,:),'b','linewidth',2)
    title(['-E2+palbo(1uM) ',proteinnamepalbo{i},])
    set(gca,'fontsize',20,'linewidth',2,'xtick',[0 1.0000 3.0000 6.0000 7.0000],...
        'xticklabel',{'0h','1d','3d','6d','7d'})
    legend({'Experiment','Prediction'},'location','NW','box','off')
    xlim([-0.2,7.2])
    ylim([0,2.5])
end
%% Proliferation
% plot of cell number of E2 control
figure
errorbar(tpro,OutE2proCounter.mean,OutE2proCounter.sem,'ro','linewidth',2);hold on
plot(tE2/24,simresultE2(end,:),'c','linewidth',2)
ylabel('Normalized cell number')
title('+E2')
set(gca,'Fontsize',20,'linewidth',2,'ytick',1:8:34,'xtick',[0 1.0000 3.0000 6.0000 7.0000],...
    'xticklabel',{'0h','1d','3d','6d','7d'})
legend({'Experiment','Simulation'},'location','NW','box','off')
ylim([0,34])
xlim([-.1,8])
grid on

% plot of cell number of +E2+ICI(500nM)
figure
errorbar(tpro,OutE2ICI500nMproCounter.mean,OutE2ICI500nMproCounter.sem,'ro','linewidth',2);hold on
plot(tE2ICI/24,simresultE2ICI(end,:),'c','linewidth',2)
ylabel('Normalized cell number')
title('+E2+ICI(500nM)')
set(gca,'Fontsize',20,'linewidth',2,'ytick',1:8:34,'xtick',[0 1.0000 3.0000 6.0000 7.0000],...
    'xticklabel',{'0h','1d','3d','6d','7d'})
legend({'Experiment','Simulation'},'location','NW','box','off')
ylim([0,34])
xlim([-.1,8])
grid on

% plot of cell number of -E2
figure
errorbar(tpro,OutE2depproCounter.mean,OutE2depproCounter.sem,'ro','linewidth',2);hold on
plot(tE2dep/24,simresultE2dep(end,:),'c','linewidth',2)
ylabel('Normalized cell number')
title('-E2')
set(gca,'Fontsize',20,'linewidth',2,'ytick',1:8:34,'xtick',[0 1.0000 3.0000 6.0000 7.0000],...
    'xticklabel',{'0h','1d','3d','6d','7d'})
legend({'Experiment','Simulation'},'location','NW','box','off')
ylim([0,34])
xlim([-.1,8])
grid on

% plot of cell number of -E2+ICI(500nM)
figure
errorbar(tpro,OutE2depICI500nMproCounter.mean,OutE2depICI500nMproCounter.sem,'ro','linewidth',2);hold on
plot(tE2depICI/24,simresultE2depICI(end,:),'b','linewidth',2)
ylabel('Normalized cell number')
title('-E2+ICI(500nM)')
set(gca,'Fontsize',20,'linewidth',2,'ytick',1:8:34,'xtick',[0 1.0000 3.0000 6.0000 7.0000],...
    'xticklabel',{'0h','1d','3d','6d','7d'})
legend({'Experiment','Prediction'},'location','NW','box','off')
ylim([0,34])
xlim([-.1,8])
grid on

% plot of cell number of +E2+palbo(1uM)
figure
errorbar(tpro,OutE2palbo1uMproCounter.mean,OutE2palbo1uMproCounter.sem,'ro','linewidth',2);hold on
plot(tE2palbo/24,simresultE2palbo(end,:),'c','linewidth',2)
ylabel('Normalized cell number')
title('+E2+palbo(1uM)')
set(gca,'Fontsize',20,'linewidth',2,'ytick',1:8:34,'xtick',[0 1.0000 3.0000 6.0000 7.0000],...
    'xticklabel',{'0h','1d','3d','6d','7d'})
legend({'Experiment','Simulation'},'location','NW','box','off')
ylim([0,34])
xlim([-.1,8])
grid on

% plot of cell number of -E2+palbo(1uM)
figure
errorbar(tpro,OutE2deppalbo1uMproCounter.mean,OutE2deppalbo1uMproCounter.sem,'ro','linewidth',2);hold on
plot(tE2deppalbo/24,simresultE2deppalbo(end,:),'b','linewidth',2)
ylabel('Normalized cell number')
title('-E2+palbo(1uM)')
set(gca,'Fontsize',20,'linewidth',2,'ytick',1:8:34,'xtick',[0 1.0000 3.0000 6.0000 7.0000],...
    'xticklabel',{'0h','1d','3d','6d','7d'})
legend({'Experiment','Prediction'},'location','NW','box','off')
ylim([0,34])
xlim([-.1,8])
grid on
