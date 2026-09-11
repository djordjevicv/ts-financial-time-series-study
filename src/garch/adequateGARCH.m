function res=adequateGARCH(a,x,m,s)
mlb=ceil(log(size(a,1)));
s2=s2sequenceGARCH(a,x,m,s);
st=sqrt(s2);

niz=a./st;

res=ljungbox2(niz,mlb);
% res=[mlb;alfa;q;z];
q=res(3);
z=res(4);
if q>z
    disp('Reject H0')
%     disp('Moguce je da postoji linearna zavisnost sa nekom od prethodnih mlb vrednosti')
    disp('The model is not adequate')
else
    disp('No significant linear dependence')
    disp('The model is adequate')
end
