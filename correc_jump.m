function [f,om,grvel]=correc_jump(grvel,tvis,phgr,om,nfa,tresh,perc,npoints)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Check dispersion curve for jumps   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
per=2*pi./om;
[ierr,trig1,ftrig]=trigger(grvel,om,nfa,tresh);
if ierr~=0
    %keep original
    grvelt = grvel;
    tvist  = tvis;
    phgrt  = phgr;
    njump = 0;
    % find all jumps
    nijmp = 0;
    for i = 1:nfa-1
        if abs(trig1(i+1)-trig1(i))>1.5
            nijmp = nijmp+1;
            ijmp(nijmp) = i;
        end
    end
    nii = 0;ii=zeros(100,1);
    for i =1:nijmp-1
        if (ijmp(i+1)-ijmp(i))<=npoints
            nii = nii +1;
            ii(nii) = i;
        end
    end
    % main loop by jumps
    if(nii~=0)
        for ki = 1:nii
            kk = ii(ki);
            grvel1 = grvelt;
            tvis1  = tvist;
            phgr1  = phgrt;
            istrt = ijmp(kk);
            ibeg = istrt+1;
            iend = ijmp(kk+1);
            ima = 0;
            for k = ibeg:iend
                dmaxt = 1e10;
                for j = 1:ici
                    if (ind(1,j)==k)
                        wor = abs(dist/(ipar(1,j))-grvel1(k-1))
                        if(wor<dmaxt)
                            ima = j
                            dmaxt = wor
                        end
                    end
                end
                grvel1(k) = dist/(ipar(1,ima))
                tvis1(k)  = ipar(2,ima)
                phgr1(k)  = ipar(6,ima)
                
                trigger(grvel1,om,nfa,tresh,trig1, ftrig1,ierr1);
                iflag = 0
                for k=istrt:iend+1
                    if(dabs(trig1(k))>=0.5)
                        iflag = 1
                    end
                end
                if(iflag==0)
                    for i =1:nfa
                        grvelt(i) = grvel1(i)
                        tvist(i)  = tvis1(i)
                        phgrt(i)  = phgr1(i)
                        njump = njump+1
                    end
                end
            end
        end
        for i=1:nfa
            grvel1(i) = grvelt(i)
            tvis1(i)  = tvist(i)
            phgr1(i)  = phgrt(i)
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % after removing possible jumps, we cut frequency range to single
        % segment with max length
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        trigger(grvel1,om,nfa,tresh,trig1, ftrig1,ierr1)
        if (ierr1~=0)
            nindx = 1
            indx(1) = 1
            for i =1:nfa
                if (abs(trig1(i))>=0.5)
                    nindx = nindx+1
                    indx(nindx) = i
                end
            end
            nindx = nindx+1
            indx(nindx) = nfa
            imax = 0
            ipos = 0
            for i =1:nindx-1
                iimax = indx(i+1)-indx(i)
                if(iimax>imax)
                    ipos = i
                    imax = iimax
                end
            end
            ist = max0(indx(ipos),1)
            ibe = min0(indx(ipos+1),nfa)
            nfout2 = ibe -ist+1
            i = ist:ibe;
            per2(i-ist+1)   = per(i);
            grvel1(i-ist+1) = grvel1(i);
            tvis1(i-ist+1)  = tvis1(i);
            phgr1(i-ist+1)  = phgr1(i);
            om1(i-ist+1)    = om(i);
            
           [ierr1,trig1,ftrig1]=trigger(grvel1,om1,nfout2,tresh);
            if(nfout2 < nfa*perc/100.0)
                ierr1 = 1;
                nfout2 = 0;
            end
        else
            nfout2 = nfa;
            per2   = per;
        end
    else
        ierr = 0;
        nfout2 = nfa;
        per2   = per;
        tvis1  = tvis;
        phgr1  = phgr;
        grvel1 = grvel;
    end
end

om=2*pi./tvis;
f=1/(2*pi)*om;