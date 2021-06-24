function [p] = Fisher(tumor,normal)

[R1,P1]=corrcoef(normal');
R1(isnan(R1))=0;
[R2,P2]=corrcoef(tumor');
R2(isnan(R2))=0;
Z1=1/2*log(abs((1+R1)./(1-R1)));
Z2=1/2*log(abs((1+R2)./(1-R2)));
dZ=abs(Z2-Z1)/(((1/(length(normal(1,:))-3)+1/(length(tumor(1,:))-3)))^0.5);
dZ(isnan(dZ))=0;
p=1-normcdf(abs(dZ));
clear  R1 P1 R2 P2 Z1 Z2 dZ
end