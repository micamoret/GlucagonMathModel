function weight_res = costfcnglumath(parameters)
global kGLUCA kg KD t sgluca_1 sgluca_2 sgluca_3 sgluca_4 sgluca_5 sgluca_6 sgluca_7 sgluca_8 sgluca_9 sgluca_10 sgluca sGLUCA dataOGTT

kGLUCA=parameters(1);
KD=parameters(2);
sgluca_1=parameters(3);
sgluca_2=parameters(4);
sgluca_3=parameters(5);
sgluca_4=parameters(6);
sgluca_5=parameters(7);
sgluca_6=parameters(8);
sgluca_7=parameters(9);
sgluca_8=parameters(10);
sgluca_9=parameters(11);
sgluca_10=parameters(12);

sgluca = [sgluca_1 sgluca_2 sgluca_3 sgluca_4 sgluca_5 sgluca_6 sgluca_7 sgluca_8 sgluca_9 sgluca_10 sgluca_10]';
sGLUCA = [t sgluca];

[a,b] = size(dataOGTT);
FSD = 0.1*ones(a,1);

[time out] = sim('glucagon_model_ins',t);


%% REGULARIZATION 

% RSS 

model_out=out(:,2);
matr=zeros(a,1);

for j=1:1
    for i=1:a
        v(i,j)=(FSD(i,j)^2)*(abs(dataOGTT(i,j)))^2;
        matr(i,j)=(dataOGTT(i,j)-model_out(i,j))/sqrt(v(i,j));
    end
end


sglucapoint=parameters(3:10);
f=find(sglucapoint>=0);
sglucapoint(f)=0;

% second derivative of sgluca
dt=diff(t);

deltaSgluca=[]; 
for i=2:10 
 deltaSgluca=[deltaSgluca;parameters(2+i)-parameters(1+i)];
end
deltaSgluca2=[]; 
for i=1:8 
 deltaSgluca2=[deltaSgluca2;deltaSgluca(1+i)-deltaSgluca(i)];
end

deltaSgluca2=deltaSgluca2./dt(1:8);

% weighted residuals
weight_res = [1*matr; 1000*kGLUCA; 1*deltaSgluca2; 10*sglucapoint'];

