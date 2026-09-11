function pz=PACFAR(rniz)

% za2=norminv(alfa/2,0,1);
za2=1.96;
T=size(rniz,1);
m=ceil(log(T));
% pac=[];
pz=0;
for p=1:m
    test=ARestimate(rniz,p);
    if test(3,1)==1
        pz=p;
    end
    
end
