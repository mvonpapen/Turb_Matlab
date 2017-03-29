clear all
gwsdir='/home/vonpapen/PhD/TEMP/GWS/';

lf=dir([gwsdir,'*04301*gws']);
n=length(lf);

for i=1:n; file{i}=[gwsdir, lf(i).name]; end

for i=1:n
    [P,S]=importpsd(file{i});
    d(i)=S.dist;
    mB(i)=S.meanB;
    B(:,i)=S.Bvec;
    vcas(:,i)=S.vcas;
    varB(i)=S.varB;
    z(i)=S.z;
    dt(i,:)=S.dt;
    f=P(:,1);
    Pa(:,i)=P(:,2);
    Pe(:,i)=P(:,3);
    Fa(:,i)=P(:,4);
    Fe(:,i)=P(:,5);
    fic(i)=mB(i)*1e-9*1.6022e-19/1.6605e-27;
end