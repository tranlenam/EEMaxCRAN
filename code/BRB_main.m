function [cbv,cup,A] = BRB_main(channel,scale,nUsers,nBS,nTx,P_bc,P_SBa,C_b)
BW = 10e6;% 10 Mhz
No = -174+10*log10(BW)-30; %dBm

P_bs = 0;
P_ac = 6.8;
P_sl = 4.3; %[mW] circuit power consumption at one BS
PAeff = 0.55; % power amplifier efficiency
P_ms = 0.1; %[W] circuit power consumption at one user
P_abh =  3.85; % [W]
P_sbh = 0.75; % [W]
P_olt = 20;

P_const = nUsers*P_ms + nBS*(P_sl+P_sbh)+P_olt+nBS*P_bs;

P_ac = P_ac + P_abh;
P_sl = P_sl + P_sbh;

P_ant = P_bc/nTx;
aeff = sqrt(P_ant)/PAeff;

channel=channel*scale;
sigma = scale*sqrt(10^((No)/10));
scalefactor = 1;

power = P_bc*scalefactor^2;
gamma_thr = (10^(0/10));
nIterations =1e6;
m=1;
A =zeros(nIterations+1,3);
iCBV = 1;
%% compute bound of variables s_b x_bk z_k u_bk
upper_bound_x = [];
lower_bound_x = [];
% binary variables x_bk and s_b
lower_bound_x = [lower_bound_x zeros(1,nUsers*nBS+nBS)];
upper_bound_x = [upper_bound_x ones(1,nUsers*nBS+nBS)];
% user rate variable z_k
lower_bound_x = [lower_bound_x log(1+gamma_thr)*ones(1,nUsers)];
for iBox=1:nUsers
    upper_bound_x = [upper_bound_x min([C_b,log(1+nBS*norm(channel(iBox,:))^2*power/sigma^2)])];
end
% power variables u_bk
lower_bound_x = [lower_bound_x zeros(1,nBS)];
upper_bound_x = [upper_bound_x nTx*sqrt(P_ant)*ones(1,nBS)];


%% check feasibility
fbhu=zeros(nBS,1);
fbhl=zeros(nBS,1);
for iBS=1:nBS
    fbhu(iBS)= upper_bound_x(nBS+(iBS-1)*nUsers+1:nBS+iBS*nUsers)*upper_bound_x(nBS+nUsers*nBS+1:nBS+nUsers*nBS+nUsers)';
    fbhl(iBS)= lower_bound_x(nBS+(iBS-1)*nUsers+1:nBS+iBS*nUsers)*lower_bound_x(nBS+nUsers*nBS+1:nBS+nUsers*nBS+nUsers)';
end
cup = sum(upper_bound_x(nBS+nUsers*nBS+1:nBS+nUsers*nBS+nUsers))/ ...
    (aeff*sum(lower_bound_x(nBS+nUsers*nBS+nUsers+1:nBS+nUsers*nBS+nUsers+nBS))/scalefactor^2+ (P_ac-P_sl)*max([1,sum(lower_bound_x(1:nBS))])+P_SBa*max([sum(fbhl.^m) sum(lower_bound_x(nBS+nUsers*nBS+1:nBS+nUsers*nBS+nUsers))]) +P_const);
cbv = sum(lower_bound_x(nBS+nUsers*nBS+1:nBS+nUsers*nBS+nUsers))/ ...
    (aeff*sum(upper_bound_x(nBS+nUsers*nBS+nUsers+1:nBS+nUsers*nBS+nUsers+nBS))/scalefactor^2+ (P_ac-P_sl)*max([1,sum(upper_bound_x(1:nBS))])+P_SBa*min([sum(fbhu.^m) C_b*nBS]) +P_const);
[~ , cobj,cbup, ~,~,~,~, ~,~,fsb] = BRB_EE_check_feasibility_mosek( channel,nUsers,nBS,nTx,sigma,scalefactor, lower_bound_x,upper_bound_x,C_b,P_ac,P_sl,P_const,P_SBa,P_ant, power, aeff,m,cbv,cup);
if fsb == 0 % feasible
    cbv = cobj;
    cup = cbup;
else
    cbv = 0;
end

%% box reduction for initial box
varupbound = upper_bound_x;
varlowbound = lower_bound_x;
etemp = zeros(nBS+nUsers*nBS+nUsers+nBS);
% alpha beta for x_bk
for iBox=1:(nBS+nUsers*nBS+nUsers+nBS)
    e_i = zeros(1,nBS+nUsers*nBS+nUsers+nBS);
    e_i(iBox) = 1;
    if iBox <= nBS+nUsers*nBS
        [ varalpha ] = BRB_bisection_alpha_discrete( varupbound,varlowbound,e_i,cbv, nUsers, nBS);
    else
        [ varalpha ] =  BRB_bisection_alpha_continuous( varupbound,varlowbound,e_i,cbv, nUsers, nBS, C_b, P_ac,P_sl,P_const,P_SBa, power, PAeff,scalefactor,m);
    end
    etemp(iBox,:) = varalpha*(varupbound(iBox)-varlowbound(iBox))*e_i;
    
end
varlowbound = varupbound - sum(etemp,1);
etemp =zeros(nBS+nUsers*nBS+nUsers+nBS);
for iBox=1:(nBS+nUsers*nBS+nUsers+nBS)
    e_i = zeros(1,nBS+nUsers*nBS+nUsers+nBS);
    e_i(iBox) = 1;
    if iBox <= nBS+nUsers*nBS
        [ varbeta ] = BRB_bisection_beta_discrete( varupbound,varlowbound,e_i,cbv, nUsers, nBS, C_b, P_ac,P_sl,P_const,P_SBa, power, PAeff,scalefactor,m);
    else
        [ varbeta ] =  BRB_bisection_beta_continuous( varupbound,varlowbound,e_i,cbv, nUsers, nBS, C_b, P_ac,P_sl,P_const,P_SBa, power, PAeff,scalefactor,m);
    end
    etemp(iBox,:) = varbeta*(varupbound(iBox)-varlowbound(iBox))*e_i;
end
varupbound = varlowbound + sum(etemp,1);

ibox = 1;
sBox = [varlowbound ;varupbound];
BOXupb(ibox) = cup;
BOXlowb(ibox) = cbv;
sboxup= cup;
sboxlow = cbv;
jstored =1;
jcand = 1;
BOXdone{jstored} = [];

%% BRB loop

for iIteration=1:nIterations  % check convergence condition
    % while iIteration <100
    %% branching
    if sum(sboxup - sboxlow)/(nUsers) >=1e-2 %check if upper and lower bound of a box converged, if true, stop processing this box
        jBox=0;
        upval =[];
        lowval =[];
        clear BOXtemp
        
        varlowboundbox1 = sBox(1,:);
        varlowboundbox2 = sBox(1,:);
        varupboundbox1 = sBox(2,:);
        varupboundbox2 = sBox(2,:);
        
        if nnz(sBox(2,1:nBS)- sBox(1,1:nBS)) % branching variable s first
            [~ ,ind] = max(sBox(2,1:nBS)- sBox(1,1:nBS));
        else
            [~ ,ind] = max(sBox(2,1:nBS+nUsers*nBS+nUsers)- sBox(1,1:nBS+nUsers*nBS+nUsers));
        end
        if ind > nBS+ nUsers*nBS
            varupboundbox1(ind) =  sBox(2,ind) - (sBox(2,ind) - sBox(1,ind))/2;
            varlowboundbox2(ind) = sBox(1,ind) + (sBox(2,ind) - sBox(1,ind))/2;
        end
        if ind <= nBS+nUsers*nBS
            
            if (sBox(1,ind) == 0) & (sBox(2,ind) == 1)
                varupboundbox1(ind) = 0;
                varlowboundbox2(ind) = 1;
            end
            
        end
        %% box reduction
        %% box 1
        etemp = zeros(1,nBS+nUsers*nBS+nUsers+nBS);
        % alpha beta for x_bk
        for iBox=1:(nBS+nUsers*nBS+nUsers+nBS)
            e_i = zeros(1,nBS+nUsers*nBS+nUsers+nBS);
            e_i(iBox) = 1;
            if iBox <= nBS+nUsers*nBS
                [ varalpha ] = BRB_bisection_alpha_discrete( varupboundbox1,varlowboundbox1,e_i,cbv, nUsers, nBS);
            else
                [ varalpha ] = BRB_bisection_alpha_continuous( varupboundbox1,varlowboundbox1,e_i,cbv, nUsers, nBS, C_b, P_ac,P_sl,P_const,P_SBa, power, PAeff,scalefactor,m);
            end
            etemp(iBox,:) = varalpha*(varupboundbox1(iBox)-varlowboundbox1(iBox))*e_i;
            
        end
        varlowboundbox1 = varupboundbox1 - sum(etemp,1);
        etemp =zeros(1,nBS+nUsers*nBS+nUsers+nBS);
        for iBox=1:(nBS+nUsers*nBS+nUsers+nBS)
            e_i = zeros(1,nBS+nUsers*nBS+nUsers+nBS);
            e_i(iBox) = 1;
            if iBox <=nBS+ nUsers*nBS
                [ varbeta ] = BRB_bisection_beta_discrete( varupboundbox1,varlowboundbox1,e_i,cbv, nUsers, nBS, C_b, P_ac,P_sl,P_const,P_SBa, power, PAeff,scalefactor,m);
            else
                [ varbeta ] = BRB_bisection_beta_continuous( varupboundbox1,varlowboundbox1,e_i,cbv, nUsers, nBS, C_b, P_ac,P_sl,P_const,P_SBa, power, PAeff,scalefactor,m);
            end
            etemp(iBox,:) = varbeta*(varupboundbox1(iBox)-varlowboundbox1(iBox))*e_i;
        end
        varupboundbox1 = varlowboundbox1 + sum(etemp,1);
        
        [ sta1 , cobj,cbup, x_bk,x_bkl,t_b,~,~,~, fsb1,time_solving1] = BRB_EE_check_feasibility_mosek( channel,nUsers,nBS,nTx,sigma,scalefactor, varlowboundbox1,varupboundbox1,C_b,P_ac,P_sl,P_const,P_SBa,P_ant, power, aeff,m,cbv,cup);
        if sta1 ==0
            if fsb1 == 0
                [cbv, idx]  = max([cbv cobj]);
                if idx ==2
                    %{ 
                    BOXcand{jcand,1} = [varlowboundbox1;varupboundbox1];
                    BOXcand{jcand,2}=cbv;
                    BOXcand{jcand,3}=iIteration;
                    %}
                    jcand = jcand+1;
                end
            end
            
            if ~nnz(varlowboundbox1(1:nBS) - varupboundbox1(1:nBS))
                varlowboundbox1(nBS+nUsers*nBS+nUsers+1:nBS+nUsers*nBS+nUsers+nBS)=[ t_b];
                varupboundbox1(nBS+1:nBS+nUsers*nBS)=[ x_bk];
                varlowboundbox1(nBS+1:nBS+nUsers*nBS)=[ x_bkl];
            end
            
            lowval = [lowval cobj];
            upval = [upval cbup];
            BOXtemp{jBox+1}=[varlowboundbox1;varupboundbox1];
            jBox=jBox+1;
            
        end
        
        
        
        %% box 2
        
        etemp = zeros(nBS+nUsers*nBS+nUsers+nBS);
        % alpha beta for x_bk
        for iBox=1:(nBS+nUsers*nBS+nUsers+nBS)
            e_i = zeros(1,nBS+nUsers*nBS+nUsers+nBS);
            e_i(iBox) = 1;
            if iBox <= nBS+nUsers*nBS
                [ varalpha ] = BRB_bisection_alpha_discrete( varupboundbox2,varlowboundbox2,e_i,cbv, nUsers, nBS);
            else
                [ varalpha ] = BRB_bisection_alpha_continuous( varupboundbox2,varlowboundbox2,e_i,cbv, nUsers, nBS, C_b, P_ac,P_sl,P_const,P_SBa, power, PAeff,scalefactor,m);
            end
            etemp(iBox,:) = varalpha*(varupboundbox2(iBox)-varlowboundbox2(iBox))*e_i;
            
        end
        varlowboundbox2 = varupboundbox2 - sum(etemp,1);
        etemp =zeros(1,nBS+nUsers*nBS+nUsers+nBS);
        for iBox=1:(nBS+nUsers*nBS+nUsers+nBS)
            e_i = zeros(1,nBS+nUsers*nBS+nUsers+nBS);
            e_i(iBox) = 1;
            if iBox <= nBS+nUsers*nBS
                [ varbeta ] = BRB_bisection_beta_discrete( varupboundbox2,varlowboundbox2,e_i,cbv, nUsers, nBS, C_b, P_ac,P_sl,P_const,P_SBa, power, PAeff,scalefactor,m);
            else
                [ varbeta ] = BRB_bisection_beta_continuous( varupboundbox2,varlowboundbox2,e_i,cbv, nUsers, nBS, C_b, P_ac,P_sl,P_const,P_SBa, power, PAeff,scalefactor,m);
            end
            etemp(iBox,:) = varbeta*(varupboundbox2(iBox)-varlowboundbox2(iBox))*e_i;
        end
        varupboundbox2 = varlowboundbox2 + sum(etemp,1);
        
        [ sta2 , cobj,cbup, x_bk,x_bkl,t_b,~, ~,~,fsb2,time_solving2] = BRB_EE_check_feasibility_mosek( channel,nUsers,nBS,nTx,sigma,scalefactor, varlowboundbox2,varupboundbox2,C_b,P_ac,P_sl,P_const,P_SBa,P_ant, power, aeff,m,cbv,cup);
        
        if sta2 ==0
            if fsb2 == 0
                [cbv, idx]  = max([cbv cobj]);
                if idx ==2
                    %{
                    BOXcand{jcand,1} = [varlowboundbox2;varupboundbox2];
                    BOXcand{jcand,2}=cbv;
                    BOXcand{jcand,3}=iIteration;
                    %}
                    jcand = jcand+1;
                end
            end
            
            if ~nnz(varlowboundbox2(1:nBS) - varupboundbox2(1:nBS))
                varlowboundbox2(nBS+nUsers*nBS+nUsers+1:nBS+nUsers*nBS+nUsers+nBS)=[ t_b];
                varupboundbox2(nBS+1:nBS+nUsers*nBS)=[ x_bk];
                varlowboundbox2(nBS+1:nBS+nUsers*nBS)=[ x_bkl];
            end
            
            lowval = [lowval cobj];
            upval = [upval cbup];
            BOXtemp{jBox+1}=[varlowboundbox2;varupboundbox2];
            jBox=jBox+1;
            
        end
        
        
        
        %% delete boxes
        if jBox > 0
            for ib = 1:length(BOXtemp)
                if cbv <= upval(ib)
                    BOX{ibox}=BOXtemp{ib};
                    BOXupb(ibox) = upval(ib);
                    BOXlowb(ibox) = lowval(ib);
                    ibox = ibox+1;
                    
                end
            end
        end
        
        cup=max(BOXupb);
        clo=min(BOXlowb);
        %{
        disp([cup cbv clo]);
        disp(iIteration);
        %}
        A(iCBV,:)=[cup cbv clo]; % save the upper bound, current best objective and lower bound for display
        iCBV = iCBV+1;
    else
        BOXdone{jstored} = sBox;
        jstored = jstored +1;
    end
    
    
    boxind = find(BOXupb >= cbv);
    if length(boxind) == 1
        sBox =  BOX{boxind};
        ibox = 1;
        if iIteration >=1000
            break;
        end
    else
        BOX =  (BOX(boxind));
        BOXupb = BOXupb(boxind);
        BOXlowb = BOXlowb(boxind);
        ibox = length(BOXupb)+1;
        [~, jbox] =max(BOXupb);
        sBox = BOX{jbox};
        BOX{jbox} = [];
        sboxup = sBox(2,1:nBS+nUsers*nBS+nUsers);
        sboxlow = sBox(1,1:nBS+nUsers*nBS+nUsers);
        BOXupb(jbox) = 0;
    end
    
    if (cup-cbv)<1e-3
        break;
    end
    
end
A(iCBV:end,:)=[];



