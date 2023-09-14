function [ myalpha ] = BRB_bisection_alpha_continuous( varupbound,varlowbound,e_i,cbv, nUsers, nBS, C_b, P_ac,P_sl,P_const,P_SBa, power, aeff,scalefactor,m)
%
a = 0;
b = 1;
fbh = zeros(nBS,1);
for iBS=1:nBS
    fbh(iBS)= varlowbound(nBS+(iBS-1)*nUsers+1:nBS+iBS*nUsers)*varlowbound(nBS+nUsers*nBS+1:nBS+nUsers*nBS+nUsers)';
end
% if (nnz((fu)<=0)== nBS) & (nnz((fbh)<=C_b)==nBS)
if (nnz((fbh)<=C_b)==nBS)
    while b-a>1e-7
        
        c = (a+b)/2;
        
        tempx = varupbound - c*(varupbound-varlowbound).*e_i;
        
        fxuz = cbv*(aeff*sum((varlowbound(nUsers*nBS+nBS+nUsers+1:nUsers*nBS+nBS+nUsers+nBS)))/scalefactor^2 + (P_ac-P_sl)*max([1,sum(varlowbound(1:nBS))])+P_SBa*sum(varlowbound(1:nBS).*(fbh.^m)) +P_const) - sum(tempx(nUsers*nBS+nBS+1:nUsers*nBS+nBS+nUsers));
        
        
        
        if  (fxuz <= 0)
            a = c;
        else b = c;
        end
        
    end
    if (c <=1e-6)
        c = 0;
    end
    if (c >=(1-1e-6))
        c = 1;
    end
    myalpha = c;
else
    myalpha = 0;
end
end

