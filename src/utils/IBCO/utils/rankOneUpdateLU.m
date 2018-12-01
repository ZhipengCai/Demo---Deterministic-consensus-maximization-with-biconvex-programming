function [AInvOut] = rankOneUpdateLU(AInvIn,u,vT)
    AInvOut = AInvIn-AInvIn*u*vT*AInvIn/(1+vT*AInvIn*u);
end