function [ mybeta ] = BRB_bisection_beta_discrete( varupbound,varlowbound,e_i,cbv, nUsers, nBS, C_b, P_ac,P_sl,P_const,P_SBa, power, aeff,scalefactor,m)
% 


tempx = varlowbound + e_i;



for iBS=1:nBS
    fbh(iBS)= tempx(nBS+(iBS-1)*nUsers+1:nBS+iBS*nUsers)*tempx(nBS+nUsers*nBS+1:nBS+nUsers*nBS+nUsers)';

end
for iBS = 1:nBS
   fsx(iBS) = (nnz(tempx(nBS+(iBS-1)*nUsers+1:nBS+iBS*nUsers) <= varupbound(iBS)) == nUsers);
end
for iBS=1:nBS

fBS(iBS) = tempx(iBS) - sum(varupbound(nBS+(iBS-1)*nUsers+1:nBS+iBS*nUsers));
end

fxuz = cbv*(aeff*sum(tempx(nUsers*nBS+nBS+nUsers+1:nUsers*nBS+nBS+nUsers+nBS))/scalefactor^2+(P_ac-P_sl)*max([1,sum(tempx(1:nBS))])+P_SBa*sum(tempx(1:nBS).*(fbh.^m)) +P_const) - sum(varupbound(nUsers*nBS+nBS+1:nUsers*nBS+nBS+nUsers));


if (nnz((fbh)<=C_b)==nBS) & (fxuz <= 0) & (nnz(fsx) == nBS) & (nnz(fBS<=0) == nBS) 
mybeta = 1;
else mybeta = 0;
end


end

