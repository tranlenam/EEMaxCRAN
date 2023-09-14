function [ mybeta ] = BRB_bisection_beta_continuous( varupbound,varlowbound,e_i,cbv, nUsers, nBS, C_b,  P_ac,P_sl,P_const,P_SBa, power, aeff,scalefactor,m)
%
a = 0;
b = 1;


while b-a>1e-7
    
    c = (a+b)/2;
    
    tempx = varlowbound + c*(varupbound-varlowbound).*e_i;
    fbh = zeros(nBS,1);
    for iBS=1:nBS
        %     fbh(iBS)= varlowbound((iBS-1)*nUsers+1:iBS*nUsers)*varlowbound(2*nUsers*nBS+1:2*nUsers*nBS+nUsers)'-C_b;
        fbh(iBS)= tempx(nBS+(iBS-1)*nUsers+1:nBS+iBS*nUsers)*tempx(nBS+nUsers*nBS+1:nBS+nUsers*nBS+nUsers)';
        %     fu(iBS) = sum(varlowbound(nUsers*nBS+(iBS-1)*nUsers+1:nUsers*nBS+iBS*nUsers))-power;
        %       fu(iBS) = sum(tempx(nBS+nUsers*nBS+nUsers+(iBS-1)*nUsers+1:nBS+nUsers*nBS+nUsers+iBS*nUsers))-varupbound(iBS)*power;
    end
    
    % for iBS=1:nBS
    % %     fbh(iBS)= tempx((iBS-1)*nUsers+1:iBS*nUsers)*tempx(2*nUsers*nBS+1:2*nUsers*nBS+nUsers)'-C_b;
    % %     fu(iBS) = sum(tempx(nUsers*nBS+(iBS-1)*nUsers+1:nUsers*nBS+iBS*nUsers))-power;
    %     fbh(iBS)= tempx(nBS+(iBS-1)*nUsers+1:nBS+iBS*nUsers)*tempx(nBS+nUsers*nBS+1:nBS+nUsers*nBS+nUsers)'-C_b;
    %     fu(iBS) = tempx(nBS+nUsers*nBS+nUsers+iBS)-power;
    % end
    % fu = sum(varlowbound(nUsers*nBS+1:2*nUsers*nBS))-tempx(nUsers*nBS+nUsers+nBS+1)+nUsers*P_ms;
    fxuz = cbv*(aeff*sum(tempx(nUsers*nBS+nBS+nUsers+1:nUsers*nBS+nBS+nUsers+nBS))...
        /scalefactor^2+(P_ac-P_sl)*max([1,sum(tempx(1:nBS))])+...
        P_SBa*sum(tempx(1:nBS).*(fbh.^m)) +P_const) ...
        - sum(varupbound(nUsers*nBS+nBS+1:nUsers*nBS+nBS+nUsers));
    
    if  (nnz((fbh)<=C_b)==nBS) & (fxuz <= 0)
        % if (nnz((fu)<=0)==nBS) & (nnz((fbh)<=C_b)==nBS) & (fxuz <= 0)
        a = c;
    else
        b = c;
    end
    
end
if (c <=1e-6)
    c = 0;
end
if (c >=(1-1e-6))
    c = 1;
end

mybeta = c;
end

