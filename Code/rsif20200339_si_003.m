% Mathematical modeling of breast cancer cells in response to endocrine therapy and Cdk4/6 inhibition
% Journal of the Royal Society Interface
% Wei He, Diane M. Demas, Isabel P. Conde, Ayesha N. Shajahan-Haq, William T. Baumann

function x_output = E2ICIpalbo_model(t,x,PAR,ICI,E2,palbo,E2dep)
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

% Model
% ER and ICI
dER = k_ER - kd_ER * ER... % translation and degradation
    - kb_E2ER * E2 * ER + kub_E2ER * E2ER...  % bound unbound with E2
    - kb_ICIER * ICI * ER + kub_ICIER * ICIER; % bound unbound with ICI

dE2ER = - kd_E2ER * E2ER... %  degradation
    + kb_E2ER * E2 * ER - kub_E2ER * E2ER; % bound unbound between E2 and ER

dICIER = kb_ICIER * ICI * ER - kub_ICIER * ICIER... % bound unbound between ICI and ER
    - kd_ICIER * ICIER;

% cyclinD1 protein
dcyclinD1 = - kd_cyclinD1 * cyclinD1... % degradation
    + k_cyclinD1 * (1 + k_cyclinD1E2ER * E2ER^p_cyclinD1E2ER_2/(p_cyclinD1E2ER_1^p_cyclinD1E2ER_2 + E2ER^p_cyclinD1E2ER_2))... %  basal and increased translation by E2ER
    - kb_cyclinD1p21 * cyclinD1 * p21 + kub_cyclinD1p21 * cyclinD1p21... % bound unbound bewteen cyclinD1 and p21
    - kb_cyclinD1palbo * cyclinD1 * palbo + kub_cyclinD1palbo * cyclinD1palbo; % bound unbound between cyclinD1 and palbo

dcyclinD1p21 = - kd_cyclinD1 * cyclinD1p21... % degradation
    + kb_cyclinD1p21 * cyclinD1 * p21 - kub_cyclinD1p21 * cyclinD1p21; % bound unbound between cyclinD1 and p21

dcyclinD1palbo = - kd_cyclinD1 * cyclinD1palbo...% degradation
    + kb_cyclinD1palbo * cyclinD1 * palbo - kub_cyclinD1palbo * cyclinD1palbo; % bound unbound between cyclinD1 and palbo

% cMyc protein
dcMyc = - kd_cMyc * cMyc... % degradation
    + k_cMyc * (1 + k_cMycE2ER * E2ER^p_cMycE2ER_2/(p_cMycE2ER_1^p_cMycE2ER_2 + E2ER^p_cMycE2ER_2)... % basal and increased translation by E2ER
    + k_cMycppRb * ppRb^p_cMycppRb_2/(p_cMycppRb_1^p_cMycppRb_2 + ppRb^p_cMycppRb_2)); % increased translation by cell cycle

% p21 protein
dp21 = - kd_p21 * p21... % transcription and degradation
    + k_p21 * (p_p21cMyc_1^p_p21cMyc_2/(p_p21cMyc_1^p_p21cMyc_2 + cMyc^p_p21cMyc_2))... % basal and inhibited transcription by cMyc
    - kb_cyclinD1p21 * cyclinD1 * p21 + kub_cyclinD1p21 * cyclinD1p21... % bound unbound bewteen cyclinD1 and p21   
    - kb_cyclinEp21 * cyclinE * p21 + kub_cyclinEp21 * cyclinEp21; % bound unbound between cyclinE and p21

% cyclinE protein
dcyclinE = k_cyclinE - kd_cyclinE * cyclinE... % translation and degradation
    - kb_cyclinEp21 * cyclinE * p21 + kub_cyclinEp21 * cyclinEp21; % bound between p21 and cyclinE
    
dcyclinEp21 = - kd_cyclinE * cyclinEp21... % degradation
    + kb_cyclinEp21 * cyclinE * p21 - kub_cyclinEp21 * cyclinEp21;% bound ubound between p21 and cyclinE

% Rb protein
dRb = - kd_Rb * Rb...  % degradation
    + k_Rb... % basal translation
    + k_RbppRb * ppRb^p_RbppRb_2/(p_RbppRb_1^p_RbppRb_2 + ppRb^p_RbppRb_2)...% increased transcription by ppRb
    - k_RbcyclinD1 * cyclinD1 * Rb^p_cyclinD1Rb_2/(p_cyclinD1Rb_1^p_cyclinD1Rb_2 + Rb^p_cyclinD1Rb_2)... % phosphorylated by cyclinD1
    + k_pRbdepho * pRb^p_pRb_2/(p_pRb_1^p_pRb_2 + pRb^p_pRb_2); % dephosphorylation of pRb

dpRb = - kd_pRb * pRb... % degradation
    + k_RbcyclinD1 * cyclinD1 * Rb^p_cyclinD1Rb_2/(p_cyclinD1Rb_1^p_cyclinD1Rb_2 + Rb^p_cyclinD1Rb_2)... % phosphorylated by cyclinD1
    - k_pRbdepho * pRb^p_pRb_2/(p_pRb_1^p_pRb_2 + pRb^p_pRb_2)... % dephosphorylation of pRb
    - k_RbcyclinE * cyclinE * pRb^p_cyclinEpRb_2/(p_cyclinEpRb_1^p_cyclinEpRb_2 + pRb^p_cyclinEpRb_2)... % phophorylated by cyclinE, related to cMyc, cMyc transcription factor of cdc25A 
    + k_ppRbdepho * ppRb^p_ppRb_2/(p_ppRb_1^p_ppRb_2 + ppRb^p_ppRb_2); % dephosphorylation of ppRb

dppRb = -kd_ppRb * ppRb... %degradation
    + k_RbcyclinE * cyclinE * pRb^p_cyclinEpRb_2/(p_cyclinEpRb_1^p_cyclinEpRb_2 + pRb^p_cyclinEpRb_2)... % phophorylated by cyclin E, related to cMyc, cMyc transcription factor of cdc25A 
    - k_ppRbdepho * ppRb^p_ppRb_2/(p_ppRb_1^p_ppRb_2 + ppRb^p_ppRb_2); % dephosphorylation of ppRb

% Cell
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