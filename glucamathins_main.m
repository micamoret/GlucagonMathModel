% glucamathins_main runs the Simulink implementation of the mathematical
% model describing glucagon kinetics as driven by insulin during an oral
% glucose tolerance test (ogtt). The model represents a modified version
% (considering insulin instead of C-peptide) of the model proposed in:
%
% Morettini M, Burattini L, Göbl C, Pacini G, Ahrén B, Tura A. Mathematical
% Model of Glucagon Kinetics for the Assessment of Insulin-Mediated
% Glucagon Inhibition During an Oral Glucose Tolerance Test. Front
% Endocrinol. 2021;12:611147. doi:10.3389/fendo.2021.611147.
%
% glucamathins_main allows to estimate individual model parameters using
% lsqnonlin to solve a nonlinear least squares problem with regularization.
%
% Following steps should be taken before running glucamathins_main:
%   1. Define in glucamathins-main the individual glucagon and insulin
%      concentrations during the ogtt and the related time samples.
%   2. Adapt the number of Sgluca values to estimate in relation to the
%      number of ogtt time samples either in glucamathins-main either in
%      costfcnglumath.
%   3. Now RUN glucamathins_main.
% 
% VERSIONS:
%   v1.0 - first version
%
% Created in Matlab R2023b.
%
% Copyright 2024 Micaela Morettini & Andrea Tura

global kGLUCA kg KD t sgluca_1 sgluca_2 sgluca_3 sgluca_4 sgluca_5 sgluca_6 sgluca_7 sgluca_8 sgluca_9 sgluca_10 sgluca sGLUCA dataOGTT 

% -------------------------------------------------------------------------
% HERE THE USER IS REQUIRED TO PROVIDE THE OGTT INFO
% definition of the individual glucagon and insulin concentrations during
% the ogtt and the related time samples.
% -------------------------------------------------------------------------
t = [0:30:300]'; % define here the ogtt time samples (min)
insulin = [...]; % define here the insulin concentrations during ogtt as a column vector (min)
glucagon = [...]; % define here the glucagon concentrations during ogtt as a column vector (min)
% -------------------------------------------------------------------------
INSb = insulin(1);
INSplasma = [t (insulin-INSb)];

GLUCAb = glucagon(1);
dataOGTT = [glucagon];
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Initial values for the parameters to be estimated. The number of
% estimated Sgluca values is equal to the number of ogtt time samples - 1.
% In this case, 11 ogtt time samples implies 10 estimated Sgluca values. 

kGLUCA_0 = 0.005;
KD_0 = 0.155;
sgluca_1_0 = 0.1;
sgluca_2_0 = 0.1;
sgluca_3_0 = 0.1;
sgluca_4_0 = 0.1;
sgluca_5_0 = 0.1;
sgluca_6_0 = 0.1;
sgluca_7_0 = 0.1;
sgluca_8_0 = 0.1;
sgluca_9_0 = 0.1;
sgluca_10_0 = 0.1;


startpar = [kGLUCA_0 KD_0 sgluca_1_0 sgluca_2_0 sgluca_3_0 sgluca_4_0 sgluca_5_0 sgluca_6_0 sgluca_7_0 sgluca_8_0 sgluca_9_0 sgluca_10_0];
sgluca_0 = [sgluca_1_0 sgluca_2_0 sgluca_3_0 sgluca_4_0 sgluca_5_0 sgluca_6_0 sgluca_7_0 sgluca_8_0 sgluca_9_0 sgluca_10_0];

% -------------------------------------

options=optimset('Display','Iter','TolFun',1E-7,'TolX',1E-6,'MaxFunEvals',3000,'MaxIter',500);

[bestpar,resnorm,residual,exitflag,output,lambda,jacobian]=lsqnonlin(@costfcnglumath,startpar,[0 0 -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf],[1 1 inf inf inf inf inf inf inf inf inf inf],options);


kGLUCA=bestpar(1);
KD=bestpar(2);
sgluca_1=bestpar(3);
sgluca_2=bestpar(4);
sgluca_3=bestpar(5);
sgluca_4=bestpar(6);
sgluca_5=bestpar(7);
sgluca_6=bestpar(8);
sgluca_7=bestpar(9);
sgluca_8=bestpar(10);
sgluca_9=bestpar(11);
sgluca_10=bestpar(12);