function [AInvOut] = rankOneUpdateLU_Ac(AInvIn,u,idx)
    AInvU = AInvIn*u;
    AInvOut = AInvIn - AInvU*AInvIn(idx,:)/(1+AInvU(idx));
end