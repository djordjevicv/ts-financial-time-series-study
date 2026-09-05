function [Rt Gt logt]=returns(r)
%r is a time series column-vector of prices
Gt=r(2:end)./r(1:end-1);
Rt=Gt-1;
logt=log(Gt);
plot(Gt,'g');
hold on
plot(Rt,'r');
hold on
plot(logt);
legend('Gt','Rt','logt')