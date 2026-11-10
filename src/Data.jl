struct Data
    insfile::String; 
    I::Int; J::Int; K::Int; L::Int; Jmax::Int; Kmax::Int; 
    Mij::Array{}; Mjk::Array{}; Mkl::Array{};
    demand::Array{}; bigM::Int; cas::Array{}; cap::Array{}; cad::Array{};  
    fc::Array{}; tcp::Array{}; tcd::Array{}; tcc::Array{}; cep::Array{}; ced::Array{}; cec::Array{}; 
    vcs::Array{}; vcp::Array{}; vcd::Array{}; ves::Array{}; vep::Array{}; ved::Array{};
    Lcapasp::Array{}; Lcapapd::Array{}; gij::Array{}; gjk::Array{}; gkl::Array{};

    function Data(insfile)
        dt = readdlm(insfile)
        I,J,K,L,Jmax,Kmax = dt[1:6]
        mijloc = findall(x->x=="Mij", dt)[1][1]
        mjkloc = findall(x->x=="Mjk", dt)[1][1]
        mklloc = findall(x->x=="Mkl", dt)[1][1]
        dloc = findall(x->x=="demand", dt)[1][1]
        casloc = findall(x->x=="cas", dt)[1][1]
        caploc = findall(x->x=="cap", dt)[1][1]
        cadloc = findall(x->x=="cad", dt)[1][1]
        fcloc = findall(x->x=="fixedcost", dt)[1][1]
        tcpl = findall(x->x=="tcp", dt)[1][1]
        tcdl = findall(x->x=="tcd", dt)[1][1]
        tccl = findall(x->x=="tcc", dt)[1][1]
        cepl = findall(x->x=="cep", dt)[1][1]
        cedl = findall(x->x=="ced", dt)[1][1]
        cecl = findall(x->x=="cec", dt)[1][1]
        vcsl = findall(x->x=="vcs", dt)[1][1]
        vcpl = findall(x->x=="vcp", dt)[1][1]
        vcdl = findall(x->x=="vcd", dt)[1][1]
        vesl = findall(x->x=="ves", dt)[1][1]
        vepl = findall(x->x=="vep", dt)[1][1]
        vedl = findall(x->x=="ved", dt)[1][1]
        Lcapaspl = findall(x->x=="Lcapasp", dt)[1][1]
        Lcapapdl = findall(x->x=="Lcapapd", dt)[1][1]
        gijl = findall(x->x=="fixedcostModesp", dt)[1][1]
        gjkl = findall(x->x=="fixedcostModepd", dt)[1][1]
        gkll = findall(x->x=="fixedcostModedc", dt)[1][1]
        #######  Assigning data 
        Mij = dt[mijloc+1:mjkloc-1,:][1:I,1:J]
        Mjk = dt[mjkloc+1:mklloc-1,:][1:J,1:K]
        Mkl = dt[mklloc+1:dloc-1,:][1:K,1:L]
        demand = dt[dloc+1:dloc+L]
        bigM = sum(demand)
        cas = dt[casloc+1:caploc-1]; cap = dt[caploc+1:cadloc-1]; cad = dt[cadloc+1:fcloc-1]
        fc = dt[fcloc+1:tcpl-1,1:2]
        tcp0 = dt[tcpl+1:tcdl-1,:]; tcd0 = dt[tcdl+1:tccl-1,:]; tcc0 = dt[tccl+1:cepl-1,:];
        cep0 = dt[cepl+1:cedl-1,:]; ced0 = dt[cedl+1:cecl-1,:]; cec0 = dt[cecl+1:vcsl-1,:];
        vcs = dt[vcsl+1:vcpl-1]; vcp = dt[vcpl+1:vcdl-1,1:2]; vcd = dt[vcdl+1:vesl-1,1:2];
        ves = dt[vesl+1:vepl-1]*10^3; vep = dt[vepl+1:vedl-1,1:2]*10^3; ved = dt[vedl+1:Lcapaspl-1,1:2]*10^3;
        Lcapasp0 = dt[Lcapaspl+1:Lcapapdl-1,:]
        Lcapapd0 = dt[Lcapapdl+1:gijl-1,:]
        gij0 = dt[gijl+1:gjkl-1,:]; gjk0 = dt[gjkl+1:gkll-1,:]; gkl0 = dt[gkll+1:end,:]
        tcp,tcd,tcc,cep,ced,cec,Lcapasp,Lcapapd,gij,gjk,gkl = [],[],[],[],[],[],[],[],[],[],[]
        ########### Adjusting the structure
        ct = 1
        for i=1:I
            tmptc = []; tmpce = []; tmpL = []; tmpg = [];
            for j=1:J
                if Mij[i,j]==1
                    push!(tmptc, [round(tcp0[ct]; digits=1)]); push!(tmpce, [cep0[ct]]);  push!(tmpL, [Lcapasp0[ct]]); push!(tmpg, [gij0[ct]]); 
                else
                    push!(tmptc, round.(tcp0[ct,1:2]; digits=1)); push!(tmpce, cep0[ct,1:2]); push!(tmpL, Lcapasp0[ct,1:2]); push!(tmpg, gij0[ct,1:2]); 
                end
                ct+=1
            end
            push!(tcp,tmptc); push!(cep,[round.(p*10^3; digits=1) for p in tmpce]);  push!(Lcapasp,tmpL); push!(gij,tmpg)
        end

        ct = 1
        for j=1:J
            tmptc = []; tmpce = []; tmpL = []; tmpg = [];
            for k=1:K
                if Mjk[j,k]==1
                    push!(tmptc, [round(tcd0[ct]; digits=1)]); push!(tmpce, [ced0[ct]]); push!(tmpL, [Lcapapd0[ct]]); push!(tmpg, [gjk0[ct]]); 
                else
                    push!(tmptc, round.(tcd0[ct,1:2]; digits=1)); push!(tmpce, ced0[ct,1:2]); push!(tmpL, Lcapapd0[ct,1:2]); push!(tmpg, gjk0[ct,1:2]);
                end
                ct+=1
            end
            push!(tcd,tmptc);   push!(ced,[round.(d*10^3; digits=1) for d in tmpce]); push!(Lcapapd,tmpL); push!(gjk,tmpg)
        end

        ct = 1
        for k=1:K
            tmptc = []; tmpce = []; tmpg = [];
            for l=1:L
                if Mkl[k,l]==1
                    push!(tmptc, [round(tcc0[ct]; digits=1)]);  push!(tmpce, [cec0[ct]]); push!(tmpg, [gkl0[ct]]); 
                else
                    push!(tmptc, round.(tcc0[ct,1:2]; digits=1)); push!(tmpce, cec0[ct,1:2]); push!(tmpg, gkl0[ct,1:2]);
                end
                ct+=1
            end
            push!(tcc,tmptc); push!(cec,[round.(c*10^3; digits=1) for c in tmpce]); push!(gkl,tmpg)
        end
        # cep2 = cep*10^3; ced2 = ced*10^3; cec2 = cec*10^3; 
        new(insfile,I,J,K,L,Jmax,Kmax,Mij,Mjk,Mkl,demand,bigM,cas,cap,cad,fc,tcp,tcd,tcc,cep,ced,cec,vcs,vcp,vcd,ves,vep,ved,Lcapasp,Lcapapd,gij,gjk,gkl)
    end
end
