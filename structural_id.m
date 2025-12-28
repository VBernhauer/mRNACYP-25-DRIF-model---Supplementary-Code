function [] = structural_id()

    syms t k_pxrmax k_r L k_pxrdeg k_mRNA k_mRNAdeg;
    syms mRNA(t) pxr(t) dmRNA(t) dpxr(t);

    % assume(pxr(0)==0)
    % assume(mRNA(0)==1)

    dpxr(t)  = k_pxrmax * exp(-k_r * t) * L * (1 - pxr(t)) - k_pxrdeg * pxr(t);
    dmRNA(t) = k_mRNA * pxr(t) + k_mRNAdeg * (1 - mRNA(t));


    dpxr_0  = subs(dpxr(t),[t,pxr(t)],[0,0])
    dmRNA_0 = subs(dmRNA(t),[t,pxr(t),mRNA(t)],[0,0,1])

    ddpxr   = diff(dpxr(t),t)
    ddmRNA  = diff(dmRNA(t),t)

    dddmRNA = diff(ddmRNA,t)


end