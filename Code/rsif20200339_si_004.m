% Mathematical modeling of breast cancer cells in response to endocrine therapy and Cdk4/6 inhibition
% Journal of the Royal Society Interface
% Wei He, Diane M. Demas, Isabel P. Conde, Ayesha N. Shajahan-Haq, William T. Baumann

% Assign values to the different measured proteins
function simresult = E2ICIpalbo_assignval(x)
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

ERt = ER + E2ER + ICIER;
cyclinD1t = cyclinD1 + cyclinD1p21 + cyclinD1palbo;
p21t = p21 + cyclinD1p21 + cyclinEp21;
Rbt = Rb + pRb + ppRb;

simresult = [];
simresult(1,:) = cMyc/cMyc(1);
simresult(2,:) = cyclinD1t/cyclinD1t(1);
simresult(3,:) = ERt/ERt(1);
simresult(4,:) = p21t/p21t(1);
simresult(5,:) = ppRb/ppRb(1);
simresult(6,:) = Rbt/Rbt(1);
simresult(7,:) = cell/cell(1); 