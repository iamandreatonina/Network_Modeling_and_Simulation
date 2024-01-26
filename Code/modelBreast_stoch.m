
%% Definition of the Breast cancer model

% we obtained the reactions from the ODEs:

% r1: 0 -> ER kER
% r2: ER -> 0 kdER
% r3: ER + E2 -> E2ER  kbE2ER
% r4: ER + ICI -> ERICI kbICIER
% r5: E2ER -> E2 + ER kubE2ER
% r6: E2ER -> 0 kdE2ER
% r7: ICIER -> ER +ICI kubICIER
% r8: ICIER -> 0 kdICIER
% r9: cycD -> 0 kdcycD
% r10: 0 -> cycD kcycD+KcycDE2ER
% r11: cycD + p21 -> cycDp21 kbcycDp21
% r12: cycD + palbo -> cycDpalbo kbcycDpalbo
% r13: cycDp21 -> cycD + p21 kubcycDp21
% r14: cycDp21 -> 0 kdcycDp21
% r15: cycDpalbo -> 0 kdcycD
% r16: cycDpalbo -> cycD + palbo kubcycDpalbo
% r17: myc -> 0 kdmyc
% r18: 0 -> myc kmyc + kE2ERmyc + kmycRB1pp
% r19: p21 -> 0 kdp21
% r20: 0 -> p21 kp21 - kmycp21
% r21: cycE -> 0 kdcycE
% r22: 0 -> cycE kcycE
% r23: cycE + p21 -> cycEp21 kbcycEp21
% r24: cycEp21 -> p21 + cycE kubcycEp21
% r25: cycEp21 -> 0 kdcycE
% r26: rb1 -> 0 kdrb1
% r27: 0 -> rb1 krb1 + krb1rb1pp 
% r28: rb1 + cycD -> rb1p krb1cycD
% r29: rb1p -> rb1 krb1pdepho
% r30: rb1p -> 0 kdrb1p
% r31: rb1p + cycE -> rb1pp krb1pcycE
% r32: rb1pp -> rb1p krb1ppdepho
% r33: rb1pp -> 0 kdrb1pp

% proliferation: cell=kpro+ kpro*kprorb1pp*rb1pp

% number of total variables: 16
% reactions: 33

% vMinus : reagents
% vPlus : products

vMinus = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0;...
0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0;...
0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0;...
0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0;...
0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0;...
0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0;...
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0;...
0 0 0 0 0 1 0 0 1 0 0 0 0 0 0 0;...
0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0;...
0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0;...
0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0;...
0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0;...
0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0;...
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0;...
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0;...
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0;...
0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0;...
0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0;...
0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0;...
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0;...
0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0;...
0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0;...
0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0;...
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1;...
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1;...
];
vPlus  = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0;...
0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0;...
1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0;...
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0;...
0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0;...
0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0;...
0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0;...
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
0 0 0 0 0 1 0 0 1 0 0 0 0 0 0 0;...
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0;...
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0;...
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0;...
0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0;...
0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0;...
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0;...
0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0;...
0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0;...
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1;...
0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0;...
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
];

% Initial state proteins: Order of proteins: ER E2ER, ICIER, cyclinD1,
% cyclinD1p21, cyclinD1palbo, cMyc, p21, cyclinEp21, Rb, pRb, ppRb
x0 = ones(16,1);
x0(2)=10^(-5);
x0(5) = 0; % position for ICIER = 0
x0(10) = 0; % position for cyclinD1palbo = 0 
x0(end) = 1;
initial=x0.';

% Experimental values for the proteins in condition of E2 deprivation
E2_dep= [1,0.750,0.875,0.828,0.209,0.188;
    1,0.982,0.945,1.044,0.860,0.921;
    1,1.141,2.03,2.371,2.342,2.692;
    1,1.0156,1.134,1.163,1.249,1.337;
    1,0.976,0.892,0.844,0.573,0.547;
    1,1.218,0.433,0.729,0.387,0.413];

t= [0, 0.1667, 1, 3, 5, 7]; % day


% we move to the stochastic reaction rates, ergo cj's, so we need the
% conversion value
Avo_Vol=6*10^(15); % Avo_Vol=N_avogadro*Cell_volume

% Parameters -> k for the deterministic reactions, c for the stochastic
% ones.
k=[0.0207; 0.1; 4266.4776;206.8202;1;0.3;1; 0.5184;1.73;11.7109;32.0864;...
5.1776;1;0.1375;0.1375;1;2.31;51.8;1.2879;5.3598;0.8058;2.4352;1.39;20.9396;0.8058;...
3.267;11.7;7;15.0134;23.6192;0.8;5.3707;9.2727];

c=[0.0207*Avo_Vol; 0.1*Avo_Vol; 4266.4776/Avo_Vol;206.8202/Avo_Vol;1*Avo_Vol;0.3/Avo_Vol;1*Avo_Vol; 0.5184*Avo_Vol;1.73*Avo_Vol;11.7109/Avo_Vol;32.0864/Avo_Vol;...
5.1776*Avo_Vol;1*Avo_Vol;0.1375/Avo_Vol;0.1375*Avo_Vol;1*Avo_Vol;2.31*Avo_Vol;51.8*Avo_Vol;1.2879*Avo_Vol;5.3598*Avo_Vol;0.8058/Avo_Vol;2.4352/Avo_Vol;1.39*Avo_Vol;20.9396*Avo_Vol;0.8058*Avo_Vol;...
3.267*Avo_Vol;11.7*Avo_Vol;7/Avo_Vol;15.0134;23.6192*Avo_Vol;0.8/Avo_Vol;5.3707;9.2727*Avo_Vol];


% chosen delta=0.1
[T, Dynamics]=simRSSA(vMinus, vPlus, k, initial, 0.1, 7, 0.1)

proliferation=kpro+(kpro*kprorb1pp*Dynamics(:,16));
plot(T, Dynamics(:,11))

figure;
plot(t,E2_dep(1,:),'-ro','linewidth',2);hold on
plot(T,Dynamics(:,11),'b','linewidth',2)
title(['-E2 ','c-Myc'])
set(gca,'fontsize',20,'linewidth',2,'xtick',[0 0.1667 1.0000 3.0000 5.0000 7.0000],...
        'xticklabel',{'','4h','1d','3d','5d','7d'})
legend({'Exp','Sim'},'location','NW','box','off')
xlim([-0.2,7])
ylim([0,35])
saveas(gca, ['ourData ','.png'])

%% RSSA IMPLEMENTATION FROM PROF. MARCHETTI - code given to us at lecture

function [T,Dynamics] = simRSSA(vMinus,vPlus,c,initialState,delta,tMax,dT)

tic % this allows to compute the simulation runtime
    
    % T and Dynamics initialization
    T = (0:dT:tMax)';
    Dynamics = nan(length(T),length(initialState)); % for each step we provide the abundance of each variable
    
    v = vPlus - vMinus;

    % variable used to limit the number of allowed simulated steps
    maxAllowedSteps = tMax*500000; % we allow at most an average of 500,000 reaction events per unit of time
    
    % pre-generation of some random numbers
    initialLenght = 1000;
    randV = rand(1,initialLenght);
    nRandVResets = 1; % I keep in memory how many times I generate the vector, so that I can compute how many random numbers have been used
    usedRandomNumbers = 0;
    
  
    % setting initial state
    i = 1;
    T(i) = 0;
    Dynamics(i,:) = initialState;
    currentTime = 0;
    currentState = initialState;
    currentStateUp = initialState + round(delta*initialState); % state upper bound
    currentStateDown = initialState - round(delta*initialState); % state lower bound
    
    % check for negative state lower bounds (in case we shift the
    % flucutation interval to have the lower bound set to zero)
    currentStateUp(currentStateDown < 0) = currentStateUp(currentStateDown < 0) + abs(currentStateDown(currentStateDown < 0));
    currentStateDown(currentStateDown < 0) = currentStateDown(currentStateDown < 0) + abs(currentStateDown(currentStateDown < 0));
    
    % computation of reaction propensities (up and down)
    aUp = zeros(size(c));
    aDown = zeros(size(c));
    for j = 1:length(c)
        aUp(j) = computeReactionPropensity(vMinus,c,currentStateUp,j);
        aDown(j) = computeReactionPropensity(vMinus,c,currentStateDown,j);
    end
    
    % sum of the upper bound of reaction propensities
    a0Up = sum(aUp);

    % simulation loop
    nSimulationSteps = 0;
    nFastAccept = 0;
    nSlowAccept = 0;
    nRejections = 0;
    nFluctuationIntUpdates = 0;
    while (currentTime < tMax && nSimulationSteps < maxAllowedSteps)
        stateConsistency = true;
        while (stateConsistency && currentTime < tMax && nSimulationSteps < maxAllowedSteps)
            u = 1;
            reactionAccepted = false;
            if (a0Up > 0) % this if clause should avoid searching for a reaction if no reactions can be selected...
                while (~reactionAccepted)
                    % extraction of three unused random numbers from randV
                    if (usedRandomNumbers + 3 > length(randV))
                        % generation of new random numbers if we reached the end of the array...
                        randV = rand(1,initialLenght);
                        usedRandomNumbers = 0;
                        nRandVResets = nRandVResets + 1;
                    end
                    r1 = randV(usedRandomNumbers+1);
                    r2 = randV(usedRandomNumbers+2);
                    r3 = randV(usedRandomNumbers+3);
                    usedRandomNumbers = usedRandomNumbers + 3;

                    % selection of a reaction candidate by using the upper bound of reaction propensities
                    mu = 1;
                    while sum(aUp(1:mu)) < r1*a0Up
                        mu = mu + 1;
                    end

                    % rejection-based strategy to accept/reject the reaction candidate
                    if (r2 <= aDown(mu)/aUp(mu))
                        reactionAccepted = true; % fast acceptance
                        nFastAccept = nFastAccept + 1;
                    else
                        % computation of the "real" propensity
                        a = computeReactionPropensity(vMinus,c,currentState,mu);
                        if (r2 <= a/aUp(mu))
                            reactionAccepted = true; % slow acceptance
                            nSlowAccept = nSlowAccept + 1;
                        else
                            nRejections = nRejections + 1; % if this is executed, the candidate reaction is rejected
                        end
                    end
                    u = u*r3; % in any case I update u to update the final computation of tau
                end

                % computation of the next tau
                tau = (1/a0Up)*log(1/u);

                % dynamics update
                currentTime = currentTime + tau;
                currentState = currentState + v(mu,:); % I apply reaction mu by means of its row of the stoichiometric matrix
                
                % saving of the current state to the dynamics timeseries if needed
                if currentTime >= T(i)+dT
                    i = i+1;
                    T(i) = currentTime;
                    Dynamics(i,:) = currentState
                end
            else
                % if nothing can be done, we just go ahead unitl tMax
                % without further updating the current state
                % PS: I save two times the state to be sure of capturing
                % the steady state condition of the model
                i = i+1;
                T(i) = currentTime;
                Dynamics(i,:) = currentState;
                currentTime = tMax;
                i = i+1;
                T(i) = currentTime;
                Dynamics(i,:) = currentState;
            end
            

            % update of the number of simulation steps
            nSimulationSteps = nSimulationSteps + 1;
            
            % check if the updated state is still consistent with its fluctuation interval
            stateConsistency = checkStateConsistency(currentState,currentStateUp,currentStateDown);
        end
        
        % update of the fluctation interval of the model state
        currentStateUp = currentState + round(delta*currentState); % state upper bound
        currentStateDown = currentState - round(delta*currentState); % state lower bound

        % check for negative state lower bounds (in case we shift the
        % flucutation interval to have the lower bound set to zero)
        currentStateUp(currentStateDown < 0) = currentStateUp(currentStateDown < 0) + abs(currentStateDown(currentStateDown < 0));
        currentStateDown(currentStateDown < 0) = currentStateDown(currentStateDown < 0) + abs(currentStateDown(currentStateDown < 0));
        
        % update of reaction propensities (up and down)
        aUp = zeros(size(c));
        aDown = zeros(size(c));
        for j = 1:length(c)
            aUp(j) = computeReactionPropensity(vMinus,c,currentStateUp,j);
            aDown(j) = computeReactionPropensity(vMinus,c,currentStateDown,j);
        end
        
        % update of the sum of the upper bound of reaction propensities
        a0Up = sum(aUp);
        
        nFluctuationIntUpdates = nFluctuationIntUpdates + 1;
    end
    
    % cut of the residual part of timeseries that remained set to NaN
    T = T(~isnan(Dynamics(:,1)));
    Dynamics = Dynamics(~isnan(Dynamics(:,1)),:);
    
    if (nSimulationSteps >= maxAllowedSteps)
        % printing of a warning message to tell to the user that the simulation has been stopped in advance
        disp(' '); % to print an empty line
        disp(['WARNING: the simulation reached the maximum allowed number of simulation steps (' num2str(maxAllowedSteps) ')!']);
        disp(' '); % to print an empty line
    end
    
    disp(['Total number of computed simulation steps: ' num2str(nSimulationSteps)]);
    disp(['Total number of used random numbers: ' num2str((nRandVResets-1)*initialLenght+usedRandomNumbers)]);
    disp(['Number of fluctuation interval updates: ' num2str(nFluctuationIntUpdates) ' (' num2str(nFluctuationIntUpdates/nSimulationSteps*100) '%)']);
    disp(['Number of fast acceptance steps: ' num2str(nFastAccept) ' (' num2str(nFastAccept/(nFastAccept+nSlowAccept+nRejections)*100) '%)']);
    disp(['Number of slow acceptance steps: ' num2str(nSlowAccept) ' (' num2str(nSlowAccept/(nFastAccept+nSlowAccept+nRejections)*100) '%)']);
    disp(['Number of rejection steps: ' num2str(nRejections) ' (' num2str(nRejections/(nFastAccept+nSlowAccept+nRejections)*100) '%)']);
    toc % this prints the simulation runtime (time elapsed from tic to toc) 
end

function a = computeReactionPropensity(vMinus,c,state,reactionIndex)
    a = c(reactionIndex);
    if (sum(vMinus(reactionIndex,:) > 0))
        for i = 1:length(state)
            % the following if clauses allow to limit the usage of the nchoosek function (needed only when vMinus(ReactionIndex,i) > 1)
            if vMinus(reactionIndex,i) == 1
                a = a*state(i);
            end
            if vMinus(reactionIndex,i) > 1
                if (state(i) >= vMinus(reactionIndex,i))
                    a = a*nchoosek(state(i),vMinus(reactionIndex,i)); % nchoosek(n,k) returns the binomial coefficient (n k)
                else
                    a = a*0; % no available reactants...
                end
            end
        end
    end
end

function test = checkStateConsistency(currentState,currentStateUp,currentStateDown)
    test = true;
    for i = 1:length(currentState)
        test = test && currentState(i) >= currentStateDown(i) && currentState(i) <= currentStateUp(i);
    end
end

