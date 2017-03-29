m0=[-0.2 1 1e-10 -1e-10 1]'; dm=zeros(5,1); m=m0;
nm=length(m);
nn=length(rhoe(1:3));
mindd=1e-4; dd=mindd+1;
lambda0=1e-20;
nite=1;
err=9e9;
smax=0;

while error>mindd
    [nite err(nite) smax]
    m=m+dm;
    n=1;
    for i=1:3
        for j=1:3
            for k=1:4
                l(n:n+nn-1)=squeeze(para(i,j,:,k,2))'/1e3;
                d(n:n+nn-1)=rhoe(1:nn)/1e3;
                lm(n:n+nn-1)=forw_rhol(m,v(i)/1e3,rhoi(j)/1e3,theta(k),l(n:n+nn-1));
                dd(n:n+nn-1)=d(n:n+nn-1)-lm(n:n+nn-1);
                G(n:n+nn-1,:)=jacrhol(m,v(i)/1e3,rhoi(j)/1e3,theta(k),l(n:n+nn-1));
                n=n+nn;
            end
        end
    end
    
    [U,S,V]=svd(G);
    smax=max(S(:));
    GG=Inf; lambda=lambda0/(10*nite);
    while isempty(find(isnan(GG),1))~=1 || isempty(find(isinf(GG),1))~=1
        lambda=lambda*2;
        if lambda==Inf, break, end
        GG=inv(G'*G+eye(nm)*lambda);
    end
    dm=GG*G'*dd(:)/100/nite;
    
    nite=nite+1;
    err(nite)=sqrt(sum((dd./d).^2)/nn);
    
    %% Check if error is greater than before
    lambda=lambda0; damp=1; m2=m;
    while err(nite)>err(nite-1)
        lambda=lambda*10;
        GG=Inf;
        while isempty(find(isnan(GG),1))~=1 || isempty(find(isinf(GG),1))~=1
            lambda=lambda*2;
            if lambda==Inf, break, end
            GG=inv(G'*G+eye(nm)*lambda);
        end
        dm=GG*G'*dd(:)/damp;
        m2=m+dm;
        n=1;
        for i=1:3
            for j=1:3
                for k=1:4
                    lm(n:n+nn-1)=forw_rhol(m2,v(i)/1e3,rhoi(j)/1e3,theta(k),l(n:n+nn-1));
                    dd2(n:n+nn-1)=d(n:n+nn-1)-lm(n:n+nn-1);
                    n=n+nn;
                end
            end
        end
        err(nite)=sqrt(sum((dd2./d).^2)/nn);
        damp=damp*10;
    end
    m=m2;
    
%     error=err(nite);
  
end
nite=nite-1;
m