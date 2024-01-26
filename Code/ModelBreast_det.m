%% Code for resolution of ODEs of:
% "Mathematical modeling of breast cancer cells in response to endocrine therapy and Cdk4/6 inhibition"
% Journal of the Royal Society Interface
% Wei He et al


% This code was obtained by collecting previous knowledge and analyzing
% different codes found in the internet, here the resources:

% - https://it.mathworks.com/matlabcentral/fileexchange/49766-simulation-of-breast-cancer
% - https://www.nature.com/articles/srep13583
% - Original Paper

% NB. once launched it takes a bit of time (at least in our computers)
% before the generation of all the needed figures


% Experimental protein values obtained from Western Blot from the original
% paper. Value capturing the protein elvels at different conditions:
% - E2 deprivation
% - ICI treatment
% - Palbociclib treatment
% - E2 deprivation + Palbociclib treatment

E2deprived = [1,0.750,0.875,0.828,0.209,0.188; 1,0.982,0.945,1.044,0.860,0.921; 1,1.141,2.03,2.371,2.342,2.692; 1,1.0156,1.134,1.163,1.249,1.337; 1,0.976,0.892,0.844,0.573,0.547; 1,1.218,0.433,0.729,0.387,0.413];

E2ICI = [1,0.645,0.482,0.350,0.183,0.158; 1,0.946,0.820,0.845,0.761,0.775; 1,0.693,0.675,0.629,0.691,0.720; 1,1.862,1.837,1.738,1.580,1.520; 1,0.952,0.836,0.663,0.603,0.727; 1,1.167,0.373,0.099,0.314,0.294];

E2depICI= [1,0.562,0.355,0.242,0.118,0.104; 1,1.082,0.862,0.789,0.653,0.678; 1,0.681,0.593,0.510,0.470,0.522; 1,1.251,1.618,1.292,1.042,1.256; 1,0.937,0.528,0.370,0.322,0.317; 1,0.963,0.300,0.119,0.048,0.057];

E2palbo= [1,0.922,0.585,0.379,0.330; 1,0.953,0.808,0.553,0.451];

E2depraivedpalbo = [1,0.613,0.346,0.101,0.089; 1,0.634,0.489,0.107,0.140];


% Cell numbers from original paper (experimental values) in the different conditions (as above +
% control condition)

E2ControlCell= [1,1.748,5.372,19.347,22.456];

E2deprivedCells = [1,1.322,3.772,8.965,9.938];

E2ICICells = [1,1.344,2.524,3.579,3.840];

E2deprivationsICICells = [1,1.269,2.325,2.945,3.103];

E2PalboCells = [1,1.569,2.558,3.362,3.421];

E2deprivationpalboCells = [1,1.144,1.460,1.230,1.427];

% protein names to make the plots

proteinname{1} = 'c-Myc';
proteinname{2} = 'total cyclinD1';
proteinname{3} = 'ER';
proteinname{4} = 'p21';
proteinname{5} = 'RB1-pp';
proteinname{6} = 'total RB1';
proteinnamepalbo{1} = 'c-Myc';
proteinnamepalbo{2} = 'RB1-pp';


time = [0, 0.1667, 1, 3, 5, 7]; % days -> to plot on the x axis of the plots
tpro = [0,1,3,6,7]; % days -> for proliferation
ConcE2 = 0.01; % E2 concentration 0.01uM
ConcICI = 0.5; % ICI concentration 0.5uM
ConcPalbo = 1; % Palbo concentration 1uM
ind_E2dep = 58; % to obtain the conditions in deprived medium
% -> it is the index to use to extract the concentration of E2 from the set of parameters
odefun = @ode23s; % we decided to use ode23s to solve the ODEs.
% We saw with Prof. Pugliese that ode23s can be very useful in stiff models


% initial value of variables -> we have 14 variables -> obtained from
% original paper
x0 = ones(14,1);
x0(3) = 0; % position for ICIER, so initial value = 0
x0(6) = 0; % position for cyclinD1palbo, so initial value= 0 
x0(end) = 1;  % position for cell value -> so initial at 1.

% Parameter used
Parameters = [0.0207;0.1;0.3;4266.4776;1;206.8202;...
    1;0.5184;0.1375;1.73;11.5734;2.6086;0.1704;32.0864;...
    1;5.1776; 1;0.33;2.31;14.1779;0.0282;6.998;37.3017;...
    5.5045;5.3598;1.2879;1.39;9.2431;2.4352;0.8058;1.39;...
    20.9396; 1;3.267;0.35;390.8695;10.7606;7;15.0134;23.6192;...
    0.8;5.3707;9.2727;0.05;14.0135;3.1884;0.6772;7;6.3544;...
    0.3503;59.8461;2.11;0.1196;266.3026;6.0056;6;37.3859;9.6676e-05;2.3e-05];

n = 1:size(E2ICI,1); % to iterate over the parameters



E2dep = 0; % concentration of E2 in the treatment E2 deprivation
E2 = ConcE2; % concentration of E2
ICI = 0; % value of ICI, if bigger than 0, means add ICI
palbo = 0; % concentration of palbo

tspan = [0,1000]; % time frame
opt = odeset('RelTol',1e-6,'AbsTol',1e-6); % optimization values, saw at lecture with Prof. Pugliese


% Model_breast at the end of the file
[t,x] = odefun(@(t,X)Model_breast(t,X,Parameters,ICI,E2,palbo,E2dep),tspan,x0,opt);

% set the steady state -> initially not implemented, we had a look at the
% original paper to better understand
x0 = x(end,:);

%% E2 control -> so E2 value is kept normal

tspan = 0:200;
x0(3) = 0;
x0(6) = 0;
x0(end) = 1;
[tE2,xE2] = odefun(@(t,X)Model_breast(t,X,Parameters,ICI,E2,palbo,E2dep),tspan,x0,opt);
% Value_total is a function to calculate the total of the proteins
resultE2 = Value_total(xE2);

%% +E2+ICI(500nM) -> ICI treatment
E2 = ConcE2; % concentration of E2
ICI = ConcICI; % concentration of ICI
palbo = 0; % concentration of palbo
[tE2ICI,xE2ICI] = odefun(@(t,X)Model_breast(t,X,Parameters,ICI,E2,palbo,E2dep),tspan,x0,opt);
resultE2ICI = Value_total(xE2ICI);

%% -E2 -> E2 deprivation
E2 = Parameters(ind_E2dep); % concentration of E2 in deprived medium
ICI = 0; % ICI concentration
palbo = 0; % concentration of palbo
E2dep = 1;
[tE2dep,xE2dep] = odefun(@(t,X)Model_breast(t,X,Parameters,ICI,E2,palbo,E2dep),tspan,x0,opt);
resultE2dep = Value_total(xE2dep);

%% -E2+ICI(500nM) -> combination of therapies
E2 = Parameters(ind_E2dep); % E2 deprived
ICI = ConcICI; % ICI concentration
palbo = 0; % concentration of palbo
[tE2depICI,xE2depICI] = odefun(@(t,X)Model_breast(t,X,Parameters,ICI,E2,palbo,E2dep),tspan,x0,opt);
resultE2depICI = Value_total(xE2depICI);

%% +E2+palbo(1uM)
E2 = ConcE2; % concentration of E2
ICI = 0; % concentration of ICI
palbo = ConcPalbo; % concentration of palbo
E2dep = 0;
[tE2palbo,xE2palbo] = odefun(@(t,X)Model_breast(t,X,Parameters,ICI,E2,palbo,E2dep),tspan,x0,opt);
resultE2palbo = Value_total(xE2palbo);

%% -E2+palbo(1uM) -> combination of therapies
E2 = Parameters(ind_E2dep); % concentration of E2 in deprived medium
ICI = 0; % concentration of ICI
E2dep = 1;
palbo = ConcPalbo; % concentration of palbo
[tE2deppalbo,xE2deppalbo] = odefun(@(t,X)Model_breast(t,X,Parameters,ICI,E2,palbo,E2dep),tspan,x0,opt);
resultE2deppalbo = Value_total(xE2deppalbo);

%% Plots
%% plots of all proteins of -E2
for i = n
    figure;
    plot(time,E2deprived(i,:),'-ro','linewidth',2); 
    % plot the experimental values
    hold on 
    plot(tE2dep/24,resultE2dep(i,:),'b','linewidth',2)
    title(['-E2 ',proteinname{i}])
    set(gca,'fontsize',20,'linewidth',2,'xtick',[0 0.1667 1.0000 3.0000 5.0000 7.0000],...
        'xticklabel',{'','4h','1d','3d','5d','7d'})
    legend({'Experiment','Simulation'},'location','NW','box','off')
    xlim([-0.2,7.2])
    ylim([0,4])
    saveas(gca, ['mE2 ',proteinname{i}, '.png'])
    i
end
%% plots of all proteins of +E2+ICI(500nM)
for i = n
    figure;
    plot(time,E2ICI(i,:),'-ro','linewidth',2);
    % plot the experimental values
    hold on 
    plot(tE2ICI/24,resultE2ICI(i,:),'b','linewidth',2)
    title(['+E2+ICI(500nM) ',proteinname{i}])
    set(gca,'fontsize',20,'linewidth',2,'xtick',[0 0.1667 1.0000 3.0000 5.0000 7.0000],...
        'xticklabel',{'','4h','1d','3d','5d','7d'})
    legend({'Experiment','Simulation'},'location','NW','box','off')
    xlim([-0.2,7.2])
    ylim([0,4])
    saveas(gca, ['pE2pICI',proteinname{i}, '.png'])
end

%% plots of all proteins of -E2+ICI(500nM)
for i = n
    figure
    plot(time,E2depICI(i,:),'-ro','linewidth',2);
    % plot the experimental values
    hold on
    plot(tE2depICI/24,resultE2depICI(i,:),'b','linewidth',2)
    title(['-E2+ICI(500nM) ',proteinname{i},])
    set(gca,'fontsize',20,'linewidth',2,'xtick',[0 0.1667 1.0000 3.0000 5.0000 7.0000],...
        'xticklabel',{'','4h','1d','3d','5d','7d'})
    legend({'Experiment','Prediction'},'location','NW','box','off')
    xlim([-0.2,7.2])
    ylim([0,2.5])
    saveas(gca, ['mE2pICI',proteinname{i}, '.png'])
end

%% plots of all proteins of +E2+palbo(1uM)
resultE2palboplot = resultE2palbo([1,5],:);
% to exract only the needed proteins,
% the one for which we have experimental values
for i = 1 : size(resultE2palboplot,1)
    figure
    plot(tpro,E2palbo(i,:),'-ro','linewidth',2);
    hold on
    plot(tE2palbo/24,resultE2palboplot(i,:),'b','linewidth',2)
    title(['+E2+palbo(1uM) ',proteinnamepalbo{i},])
    set(gca,'fontsize',20,'linewidth',2,'xtick',[0 1.0000 3.0000 6.0000 7.0000],...
        'xticklabel',{'0h','1d','3d','6d','7d'})
    legend({'Experiment','Simulation'},'location','NW','box','off')
    xlim([-0.2,7.2])
    ylim([0,2.5])
    saveas(gca, ['pE2palbo',proteinname{i}, '.png'])
end

%% plots of all proteins of -E2+palbo(1uM)
resultE2deppalboplot = resultE2deppalbo([1,5],:);
for i = 1 : size(resultE2deppalboplot,1)
    figure
    plot(tpro,E2depraivedpalbo(i,:),'-ro','linewidth',2);hold on
    plot(tE2deppalbo/24,resultE2deppalboplot(i,:),'b','linewidth',2)
    title(['-E2+palbo(1uM) ',proteinnamepalbo{i},])
    set(gca,'fontsize',20,'linewidth',2,'xtick',[0 1.0000 3.0000 6.0000 7.0000],...
        'xticklabel',{'0h','1d','3d','6d','7d'})
    legend({'Experiment','Prediction'},'location','NW','box','off')
    xlim([-0.2,7.2])
    ylim([0,2.5])
    saveas(gca, ['mE2palbo',proteinname{i}, '.png'])
end

%% Proliferation plots
%% plot of cell number of +E2+ICI(500nM)
figure
plot(tpro,E2ICICells,'ro','linewidth',2);
hold on
plot(tE2ICI/24,resultE2ICI(end,:),'b','linewidth',2)
ylabel('Normalized cell number')
title('+E2+ICI(500nM)')
set(gca,'Fontsize',20,'linewidth',2,'ytick',1:8:34,'xtick',[0 1.0000 3.0000 6.0000 7.0000],...
    'xticklabel',{'0h','1d','3d','6d','7d'})
legend({'Experiment','Simulation'},'location','NW','box','off')
ylim([0,34])
xlim([-.1,8])
grid on
saveas(gca, 'Cells_pE2_ici.png')

%% plot of cell number of -E2
figure
plot(tpro,E2deprivedCells,'ro','linewidth',2);hold on
plot(tE2dep/24,resultE2dep(end,:),'b','linewidth',2)
ylabel('Normalized cell number')
title('-E2')
set(gca,'Fontsize',20,'linewidth',2,'ytick',1:8:34,'xtick',[0 1.0000 3.0000 6.0000 7.0000],...
    'xticklabel',{'0h','1d','3d','6d','7d'})
legend({'Experiment','Simulation'},'location','NW','box','off')
ylim([0,34])
xlim([-.1,8])
grid on
saveas(gca, 'Cells_mE2.png')

%% plot of cell number of -E2+ICI(500nM)
figure
plot(tpro,E2deprivationsICICells,'ro','linewidth',2);hold on
plot(tE2depICI/24,resultE2depICI(end,:),'b','linewidth',2)
ylabel('Normalized cell number')
title('-E2+ICI(500nM)')
set(gca,'Fontsize',20,'linewidth',2,'ytick',1:8:34,'xtick',[0 1.0000 3.0000 6.0000 7.0000],...
    'xticklabel',{'0h','1d','3d','6d','7d'})
legend({'Experiment','Prediction'},'location','NW','box','off')
ylim([0,34])
xlim([-.1,8])
grid on
saveas(gca, 'Cells_mE2_ici.png')

% plot of cell number of +E2+palbo(1uM)
figure
plot(tpro,E2PalboCells,'ro','linewidth',2);hold on
plot(tE2palbo/24,resultE2palbo(end,:),'b','linewidth',2)
ylabel('Normalized cell number')
title('+E2+palbo(1uM)')
set(gca,'Fontsize',20,'linewidth',2,'ytick',1:8:34,'xtick',[0 1.0000 3.0000 6.0000 7.0000],...
    'xticklabel',{'0h','1d','3d','6d','7d'})
legend({'Experiment','Simulation'},'location','NW','box','off')
ylim([0,34])
xlim([-.1,8])
grid on
saveas(gca, 'Cells_pE2_palbo.png')

%% plot of cell number of -E2+palbo(1uM)
figure
plot(tpro,E2deprivationpalboCells,'ro','linewidth',2);hold on
plot(tE2deppalbo/24,resultE2deppalbo(end,:),'b','linewidth',2)
ylabel('Normalized cell number')
title('-E2+palbo(1uM)')
set(gca,'Fontsize',20,'linewidth',2,'ytick',1:8:34,'xtick',[0 1.0000 3.0000 6.0000 7.0000],...
    'xticklabel',{'0h','1d','3d','6d','7d'})
legend({'Experiment','Prediction'},'location','NW','box','off')
ylim([0,34])
xlim([-.1,8])
grid on
saveas(gca, 'Cells_mE2_palbo.png')


%% Defition of the function for the ode solver -> description of the ODEs plus parameters

function x_output = Model_breast(t,x,PAR,ICI,E2,palbo,E2dep)
% variables
ER = x(1);
E2ER = x(2);
ICIER = x(3);
cyclinD1 = x(4);
cyclinD1p21 = x(5);
cyclinD1palbo = x(6);
cMyc = x(7);
p21 = x(8);
cyclinE = x(9);
cyclinEp21 = x(10);
Rb = x(11);
pRb = x(12);
ppRb = x(13);
cell = x(14);

% to understand if we are in a condition of E2 deprivation -> had to look
% at the original paper to understand
if E2dep && t/24 >= 3
   E2 = PAR(end);
end

% disp(t)

% parameter
% E2, ER, and ICI
k_ER = PAR(1); % translation rate of ER, t=0.19h --> 1/t=5.26
kd_ER = PAR(2); % degradation rate of ER, half-life: t=5h --> -1/t*log(1/2) no ligand
kd_E2ER = PAR(3); % degradation rate of E2ER, half-life: t=3h --> -1/tlog(1/2) E2 ligand
kb_E2ER = PAR(4); % bound rate between E2 and ER
kub_E2ER = PAR(5); % unbound rate between E2 and ER
kb_ICIER = PAR(6); % bound rate between ICI and ER
kub_ICIER = PAR(7); % unbound rate between ICI and ER
kd_ICIER = PAR(8); % degradation rate of ICIER

% cyclinD1 protein
k_cyclinD1 = PAR(9); % translation rate of cyclinD1, t=0.0095h --> 1/t=105.26
kd_cyclinD1 = PAR(10); % degradation rate of cyclinD1, half-life: t=0.4h --> -1/t*log(1/2)=1.7392
k_cyclinD1E2ER = PAR(11) ; % translation rate of cyclinD1 by E2ER
p_cyclinD1E2ER_1 = PAR(12); % parameter of cyclinD1 translation by E2ER
p_cyclinD1E2ER_2 = PAR(13); % parameter of cyclinD1 translation by E2ER
kb_cyclinD1p21 = PAR(14); % bound rate beteween cyclinD1 and p21
kub_cyclinD1p21 = PAR(15); % ubound rate between cyclinD1 and p21
kb_cyclinD1palbo = PAR(16); % bound rate between cyclinD1 and palbo
kub_cyclinD1palbo = PAR(17); % unbound rate between cyclinD1 and palbo

% cMyc protein
k_cMyc = PAR(18); % translation rate of cMyc, t=0.0137h --> 1/t=72.9927
kd_cMyc = PAR(19); % degradation rate of cMyc, half-life: t=0.3h --> -1/t*log(1/2)=2.3105
k_cMycE2ER = PAR(20); % translation rate of cMyc by E2ER
p_cMycE2ER_1 = PAR(21); % parameter of cMyc translation by E2ER
p_cMycE2ER_2 = PAR(22); % parameter of cMyc translation by E2ER
k_cMycppRb = PAR(23); % translation rate of cMyc by ppRb
p_cMycppRb_1 = PAR(24); % parameter of cMyc translation by ppRb
p_cMycppRb_2 = PAR(25); % parameter of cMyc translation by ppRb

% p21 protein
k_p21 = PAR(26); % translation rate of p21, t=0.0051h --> 1/t=196.0784
kd_p21 = PAR(27); % degradation rate of p21, half-life: t=0.33-1h --> -1/t*log(1/2)=1.3863
p_p21cMyc_1 = PAR(28); % parameter of p21 translation inhibited by cMyc
p_p21cMyc_2 = PAR(29); % parameter of p21 translation inhibitied by cMyc

% cyclin E protein
k_cyclinE = PAR(30); % translation rate of cyclinE, t=0.032h --> 1/t=31.25
kd_cyclinE = PAR(31); % degradation rate of cyclinE, half-life=0.5h --> -1/t*log(1/2)=1.3863
kb_cyclinEp21 = PAR(32); % bound rate with p21
kub_cyclinEp21 = PAR(33); % unbound rate with p21

% Rb protein
k_Rb = PAR(34); % translation rate of Rb, t=0.029h --> 1/t=34.48
kd_Rb = PAR(35); % degradation rate of Rb, half-life=2h --> -1/t*log(1/2)=0.3466
k_RbppRb = PAR(36); % rate of Rb transcription increased by ppRb
p_RbppRb_1 = PAR(37); % parameter of Rb transcription increased by ppRb
p_RbppRb_2 = PAR(38); % parameter of Rb transcription increased by ppRb
k_RbcyclinD1 = PAR(39); % rate of Rb phosphorylation by cyclinD1
k_pRbdepho = PAR(40); % dephosphorylation rate of pRb
kd_pRb = PAR(41); % degradation rate of pRb, half-life=2h --> -1/t*log(1/2)=0.3466
k_RbcyclinE = PAR(42); % rate cyclinE phosphorylates of Rb
k_ppRbdepho = PAR(43); % dephosphorylation rate of ppRb
kd_ppRb = PAR(44); % degradation rate of ppRb, half-life=8h --> -1/t*log(1/2)=0.0866
p_cyclinD1Rb_1 = PAR(45); % parameter of Rb phosphorylation by cyclinD1
p_cyclinD1Rb_2 = PAR(46); % parameter of Rb phosphorylation by cyclinD1
p_pRb_1 = PAR(47); % parameter of pRb dephosphorylation
p_pRb_2 = PAR(48); % parameter of pRb dephosphorylation
p_cyclinEpRb_1 = PAR(49); % parameter of pRb phosphorylation by cyclinE
p_cyclinEpRb_2 = PAR(50); % parameter of pRb phosphorylation by cyclinE
p_ppRb_1 = PAR(51); % parameter of ppRb dephosphorylation
p_ppRb_2 = PAR(52); % parameter of ppRb dephosphorylation

% proliferation
k_pro = PAR(53); % basal proliferation rate
k_proppRb = PAR(54); % proliferation increased by ppRb
p_proppRb_1 = PAR(55); % parameter of proliferation increased by ppRb
p_proppRb_2 = PAR(56); % parameter of proliferation increased by ppRb
k_carrying = PAR(57); % carrying capactity

% Model -> rewrite the ODE as saw in the paper
% ER and ICI
dER = k_ER - kd_ER * ER...
    - kb_E2ER * E2 * ER + kub_E2ER * E2ER...  
    - kb_ICIER * ICI * ER + kub_ICIER * ICIER;

dE2ER = - kd_E2ER * E2ER... 
    + kb_E2ER * E2 * ER - kub_E2ER * E2ER; 

dICIER = kb_ICIER * ICI * ER - kub_ICIER * ICIER... 
    - kd_ICIER * ICIER;

% cyclinD1 protein
dcyclinD1 = - kd_cyclinD1 * cyclinD1... 
    + k_cyclinD1 * (1 + k_cyclinD1E2ER * E2ER^p_cyclinD1E2ER_2/(p_cyclinD1E2ER_1^p_cyclinD1E2ER_2 + E2ER^p_cyclinD1E2ER_2))... 
    - kb_cyclinD1p21 * cyclinD1 * p21 + kub_cyclinD1p21 * cyclinD1p21... 
    - kb_cyclinD1palbo * cyclinD1 * palbo + kub_cyclinD1palbo * cyclinD1palbo; 

dcyclinD1p21 = - kd_cyclinD1 * cyclinD1p21... 
    + kb_cyclinD1p21 * cyclinD1 * p21 - kub_cyclinD1p21 * cyclinD1p21; 

dcyclinD1palbo = - kd_cyclinD1 * cyclinD1palbo...
    + kb_cyclinD1palbo * cyclinD1 * palbo - kub_cyclinD1palbo * cyclinD1palbo; 

% cMyc protein
dcMyc = - kd_cMyc * cMyc... 
    + k_cMyc * (1 + k_cMycE2ER * E2ER^p_cMycE2ER_2/(p_cMycE2ER_1^p_cMycE2ER_2 + E2ER^p_cMycE2ER_2)... 
    + k_cMycppRb * ppRb^p_cMycppRb_2/(p_cMycppRb_1^p_cMycppRb_2 + ppRb^p_cMycppRb_2)); 

% p21 protein
dp21 = - kd_p21 * p21... 
    + k_p21 * (p_p21cMyc_1^p_p21cMyc_2/(p_p21cMyc_1^p_p21cMyc_2 + cMyc^p_p21cMyc_2))... 
    - kb_cyclinD1p21 * cyclinD1 * p21 + kub_cyclinD1p21 * cyclinD1p21... 
    - kb_cyclinEp21 * cyclinE * p21 + kub_cyclinEp21 * cyclinEp21;

% cyclinE protein
dcyclinE = k_cyclinE - kd_cyclinE * cyclinE... 
    - kb_cyclinEp21 * cyclinE * p21 + kub_cyclinEp21 * cyclinEp21; 
    
dcyclinEp21 = - kd_cyclinE * cyclinEp21... 
    + kb_cyclinEp21 * cyclinE * p21 - kub_cyclinEp21 * cyclinEp21;

% Rb protein
dRb = - kd_Rb * Rb...  
    + k_Rb... 
    + k_RbppRb * ppRb^p_RbppRb_2/(p_RbppRb_1^p_RbppRb_2 + ppRb^p_RbppRb_2)...
    - k_RbcyclinD1 * cyclinD1 * Rb^p_cyclinD1Rb_2/(p_cyclinD1Rb_1^p_cyclinD1Rb_2 + Rb^p_cyclinD1Rb_2)... 
    + k_pRbdepho * pRb^p_pRb_2/(p_pRb_1^p_pRb_2 + pRb^p_pRb_2);

dpRb = - kd_pRb * pRb... 
    + k_RbcyclinD1 * cyclinD1 * Rb^p_cyclinD1Rb_2/(p_cyclinD1Rb_1^p_cyclinD1Rb_2 + Rb^p_cyclinD1Rb_2)... 
    - k_pRbdepho * pRb^p_pRb_2/(p_pRb_1^p_pRb_2 + pRb^p_pRb_2)... 
    - k_RbcyclinE * cyclinE * pRb^p_cyclinEpRb_2/(p_cyclinEpRb_1^p_cyclinEpRb_2 + pRb^p_cyclinEpRb_2)... 
    + k_ppRbdepho * ppRb^p_ppRb_2/(p_ppRb_1^p_ppRb_2 + ppRb^p_ppRb_2); 

dppRb = -kd_ppRb * ppRb... 
    + k_RbcyclinE * cyclinE * pRb^p_cyclinEpRb_2/(p_cyclinEpRb_1^p_cyclinEpRb_2 + pRb^p_cyclinEpRb_2)... 
    - k_ppRbdepho * ppRb^p_ppRb_2/(p_ppRb_1^p_ppRb_2 + ppRb^p_ppRb_2); 

% Cell -> formula from the original paper
dcell = k_pro/1e3 * (1 + k_proppRb * ppRb^p_proppRb_2/(p_proppRb_1^p_proppRb_2 + ppRb^p_proppRb_2)) * cell * (1 - cell/k_carrying);

x_output(1) = dER;
x_output(2) = dE2ER;
x_output(3) = dICIER;
x_output(4) = dcyclinD1;
x_output(5) = dcyclinD1p21; 
x_output(6) = dcyclinD1palbo;
x_output(7) = dcMyc;
x_output(8) = dp21;
x_output(9) = dcyclinE;
x_output(10) = dcyclinEp21;
x_output(11) = dRb;
x_output(12) = dpRb;
x_output(13) = dppRb;
x_output(14) = dcell;
x_output = x_output(:);
end

%% Assign total values -> total protein = protein free + protein bounded
function result = Value_total(x)
ER = x(:,1); 
E2ER = x(:,2);
ICIER = x(:,3);
cyclinD1 = x(:,4);
cyclinD1p21 = x(:,5);
cyclinD1palbo = x(:,6);
cMyc = x(:,7);
p21 = x(:,8);
cyclinEp21 = x(:,10);
Rb = x(:,11);
pRb = x(:,12);
ppRb = x(:,13);
cell = x(:,14);

% total protein = protein free + protein bounded
ERt = ER + E2ER + ICIER;
cyclinD1t = cyclinD1 + cyclinD1p21 + cyclinD1palbo;
p21t = p21 + cyclinD1p21 + cyclinEp21;
Rbt = Rb + pRb + ppRb;

% normalization
result = [];
result(1,:) = cMyc/cMyc(1);
result(2,:) = cyclinD1t/cyclinD1t(1);
result(3,:) = ERt/ERt(1);
result(4,:) = p21t/p21t(1);
result(5,:) = ppRb/ppRb(1);
result(6,:) = Rbt/Rbt(1);
result(7,:) = cell/cell(1); 
end