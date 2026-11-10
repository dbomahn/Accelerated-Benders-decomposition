function NDfilter(Pobj)
    copyobj = Dict();
    for i=1:length(Pobj)
        copyobj[i] = Pobj[i]
    end
    for i=1:length(Pobj)-1
        for j=i+1:length(Pobj)
            if all(round.(Pobj[i];digits=4) .>= round.(Pobj[j];digits=4)) == true #dominated by PF[j]
                copyobj[i]=nothing; break
            elseif all(round.(Pobj[j];digits=4) .>= round.(Pobj[i];digits=4)) == true
                copyobj[j]=nothing; 
            end
        end
    end
    finalobj = filter!(a->a!=nothing, collect(values(copyobj)))
    return sort!(finalobj)
end
function dominated(x,P)
    st = false
    for k=1:length(P)
        if all( x .>= P[k])# && any(x > P[k])
            st=true; break
        else
            continue
        end
    end
    return st
end

function Classic_lazy_callback(cb_data)
    yb = callback_value.(cb_data, mas.y);    ubij = callback_value.(cb_data, mas.uij);    ubjk = callback_value.(cb_data, mas.ujk)
    ubkl = callback_value.(cb_data, mas.ukl);    θb = callback_value(cb_data, mas.θ)
    subp = solve_subp(dsp,yb,ubij,ubjk,ubkl);    
    if subp.res == :OptimalityCut
        # @info "Optimality cut found"
        if θb >= round(subp.obj; digits=num)
            return
        else
            cut = @build_constraint( mas.θ ≥ dot(sm.demand,subp.α5)-dot(sm.cas,subp.α6)-
                sum(sm.cap[j]*subp.α7[j,t]*mas.y[j,t] for j=1:sm.J for t=1:2) -
                sum(sm.cad[k]*subp.α8[k,t]*mas.y[k+sm.J,t] for k=1:sm.K for t=1:2) +
                sum(sm.Lcapasp[i][j][m]*subp.α9[i,j,m]*mas.uij[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]) +
                sum(sm.Lcapapd[j][k][m]*subp.α10[j,k,m]*mas.ujk[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]) -
                sum(sm.bigM*mas.uij[i,j,m]*subp.α12[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]) -
                sum(sm.bigM*mas.ujk[j,k,m]*subp.α13[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]) -
                sum(sm.bigM*mas.ukl[k,l,m]*subp.α14[k,l,m] for k=1:sm.K for l=1:sm.L for m=1:sm.Mkl[k,l])
                )
            # cut = @build_constraint( mas.θ ≥ dot(sm.demand,subp.α5)-dot(sm.cas,subp.α6)-
            #     sum(sm.cap[j]*subp.α7[j,t]*mas.y[j,t] for j=1:sm.J for t=1:2) -
            #     sum(sm.cad[k]*subp.α8[k,t]*mas.y[k+sm.J,t] for k=1:sm.K for t=1:2) +
            #     sum(sm.Lcapasp[i][j][2]*subp.α9[i,j,2]*mas.uij[i,j,2] for i=1:sm.I for j=1:sm.J if sm.Mij[i,j]==2 ) +
            #     sum(sm.Lcapapd[j][k][2]*subp.α10[j,k,2]*mas.ujk[j,k,2] for j=1:sm.J for k=1:sm.K if sm.Mjk[j,k]==2) 
            #     # sum(sm.Lcapasp[i][j][m]*subp.α9[i,j,m]*mas.uij[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]) +
            #     # sum(sm.Lcapapd[j][k][m]*subp.α10[j,k,m]*mas.ujk[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]) 
            #     - sum(sm.bigM*mas.uij[i,j,m]*subp.α12[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]) -
            #     sum(sm.bigM*mas.ujk[j,k,m]*subp.α13[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]) -
            #     sum(sm.bigM*mas.ukl[k,l,m]*subp.α14[k,l,m] for k=1:sm.K for l=1:sm.L for m=1:sm.Mkl[k,l])
            #     )
            MOI.submit(mas.m, MOI.LazyConstraint(cb_data), cut); push!(tmpno,1)
            # push!(ocuts, (obj=subp.obj,α5=subp.α5,α6=subp.α6,α7=subp.α7,α8=subp.α8,α9=subp.α9,α10=subp.α10,α12=subp.α12,α13=subp.α13,α14=subp.α14)) # push!(tmpno,1) 
        end
    # else
    elseif subp.res == :FeasibilityCut
        # @info "Feasibility cut found"
        # cut = @build_constraint( 0 ≥ dot(sm.demand,subp.α5)-dot(sm.cas,subp.α6)-
        #     sum(sm.cap[j]*subp.α7[j,t]*mas.y[j,t] for j=1:sm.J for t=1:2) -
        #     sum(sm.cad[k]*subp.α8[k,t]*mas.y[k+sm.J,t] for k=1:sm.K for t=1:2) +
        #     # sum(sm.Lcapasp[i][j][m]*subp.α9[i,j,m]*mas.uij[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]) +
        #     # sum(sm.Lcapapd[j][k][m]*subp.α10[j,k,m]*mas.ujk[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]) 
        #     sum(sm.Lcapasp[i][j][2]*subp.α9[i,j,2]*mas.uij[i,j,2] for i=1:sm.I for j=1:sm.J if sm.Mij[i,j]==2 ) +
        #     sum(sm.Lcapapd[j][k][2]*subp.α10[j,k,2]*mas.ujk[j,k,2] for j=1:sm.J for k=1:sm.K if sm.Mjk[j,k]==2) 
        #     - sum(sm.bigM*mas.uij[i,j,m]*subp.α12[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]) -
        #     sum(sm.bigM*mas.ujk[j,k,m]*subp.α13[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]) -
        #     sum(sm.bigM*mas.ukl[k,l,m]*subp.α14[k,l,m] for k=1:sm.K for l=1:sm.L for m=1:sm.Mkl[k,l])
        #     )
        cut = @build_constraint( 0 ≥ dot(sm.demand,subp.α5)-dot(sm.cas,subp.α6)-
            sum(sm.cap[j]*subp.α7[j,t]*mas.y[j,t] for j=1:sm.J for t=1:2) -
            sum(sm.cad[k]*subp.α8[k,t]*mas.y[k+sm.J,t] for k=1:sm.K for t=1:2) +
            sum(sm.Lcapasp[i][j][m]*subp.α9[i,j,m]*mas.uij[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]) +
            sum(sm.Lcapapd[j][k][m]*subp.α10[j,k,m]*mas.ujk[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]) -
            sum(sm.bigM*mas.uij[i,j,m]*subp.α12[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]) -
            sum(sm.bigM*mas.ujk[j,k,m]*subp.α13[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]) -
            sum(sm.bigM*mas.ukl[k,l,m]*subp.α14[k,l,m] for k=1:sm.K for l=1:sm.L for m=1:sm.Mkl[k,l])
            )
        MOI.submit(mas.m, MOI.LazyConstraint(cb_data), cut); push!(tmpnf,1)
        # push!(fcuts, (α5=subp.α5,α6=subp.α6,α7=subp.α7,α8=subp.α8,α9=subp.α9,α10=subp.α10,α12=subp.α12,α13=subp.α13,α14=subp.α14))            
        ####################    When dual values are used       #######################
        # cut = @build_constraint( 0 ≥  
        #     -sum(sm.cap[j]*subp.con4[j,t]*mas.y[j,t] for j=1:sm.J for t=1:2) -
        #     sum(sm.cad[k]*subp.con5[k,t]*mas.y[k+sm.J,t] for k=1:sm.K for t=1:2) +
        #     sum((sm.Lcapasp[i][j][m]-sm.bigM)*subp.con1[i,j,m]*mas.uij[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]) +
        #     sum((sm.Lcapapd[j][k][m]-sm.bigM)*subp.con2[j,k,m]*mas.ujk[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]) -            
        #     sum(sm.bigM*mas.ukl[k,l,m]*subp.con3[k,l,m] for k=1:sm.K for l=1:sm.L for m=1:sm.Mkl[k,l])
        #     )
        
        # push!(fcuts, (α5=subp.α5,α6=subp.α6, con1 = subp.con1, con2 = subp.con2, con3 = subp.con3, con4 = subp.con4, con5 = subp.con5))
        # Relaxing the cut
        # cut = @build_constraint( 0 ≥ dot(sm.demand,floor.(subp.α5))-dot(sm.cas,floor.(subp.α6))-
            # sum(sm.cap[j]*floor.(subp.α7[j,t])*mas.y[j,t] for j=1:sm.J for t=1:2) -
            # sum(sm.cad[k]*floor.(subp.α8[k,t])*mas.y[k+sm.J,t] for k=1:sm.K for t=1:2) +
            # sum(sm.Lcapasp[i][j][2]*floor.(subp.α9[i,j,2])*mas.uij[i,j,2] for i=1:sm.I for j=1:sm.J if sm.Mij[i,j]==2 ) +
            # sum(sm.Lcapapd[j][k][2]*floor.(subp.α10[j,k,2])*mas.ujk[j,k,2] for j=1:sm.J for k=1:sm.K if sm.Mjk[j,k]==2) 
            # - sum(sm.bigM*mas.uij[i,j,m]*floor.(subp.α12[i,j,m]) for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]) -
            # sum(sm.bigM*mas.ujk[j,k,m]*floor.(subp.α13[j,k,m]) for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]) -
            # sum(sm.bigM*mas.ukl[k,l,m]*floor.(subp.α14[k,l,m]) for k=1:sm.K for l=1:sm.L for m=1:sm.Mkl[k,l]) 
            # )
    end
    push!(tmpit,1)
    return 
end

function Classic_lazy_callback_trunc(cb_data)
    yb = callback_value.(cb_data, mas.y);    ubij = callback_value.(cb_data, mas.uij);    ubjk = callback_value.(cb_data, mas.ujk)
    ubkl = callback_value.(cb_data, mas.ukl);    θb = callback_value(cb_data, mas.θ)
    subp = solve_subp(dsp,yb,ubij,ubjk,ubkl);    
    if subp.res == :OptimalityCut
        if θb >= trunc(subp.obj; digits=num) 
            return
        else
            cut = @build_constraint( mas.θ ≥ dot(sm.demand,subp.α5)-dot(sm.cas,subp.α6)-
                sum(sm.cap[j]*subp.α7[j,t]*mas.y[j,t] for j=1:sm.J for t=1:2) -
                sum(sm.cad[k]*subp.α8[k,t]*mas.y[k+sm.J,t] for k=1:sm.K for t=1:2) +
                sum(sm.Lcapasp[i][j][m]*subp.α9[i,j,m]*mas.uij[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]) +
                sum(sm.Lcapapd[j][k][m]*subp.α10[j,k,m]*mas.ujk[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]) -
                sum(sm.bigM*mas.uij[i,j,m]*subp.α12[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]) -
                sum(sm.bigM*mas.ujk[j,k,m]*subp.α13[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]) -
                sum(sm.bigM*mas.ukl[k,l,m]*subp.α14[k,l,m] for k=1:sm.K for l=1:sm.L for m=1:sm.Mkl[k,l])
                )
            MOI.submit(mas.m, MOI.LazyConstraint(cb_data), cut); push!(tmpno,1)
            # push!(ocuts, (obj=subp.obj,α5=subp.α5,α6=subp.α6,α7=subp.α7,α8=subp.α8,α9=subp.α9,α10=subp.α10,α12=subp.α12,α13=subp.α13,α14=subp.α14)) # push!(tmpno,1) 
        end
    # else
    elseif subp.res == :FeasibilityCut

        cut = @build_constraint( 0 ≥ dot(sm.demand,subp.α5)-dot(sm.cas,subp.α6)-
            sum(sm.cap[j]*subp.α7[j,t]*mas.y[j,t] for j=1:sm.J for t=1:2) -
            sum(sm.cad[k]*subp.α8[k,t]*mas.y[k+sm.J,t] for k=1:sm.K for t=1:2) +
            sum(sm.Lcapasp[i][j][m]*subp.α9[i,j,m]*mas.uij[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]) +
            sum(sm.Lcapapd[j][k][m]*subp.α10[j,k,m]*mas.ujk[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]) -
            sum(sm.bigM*mas.uij[i,j,m]*subp.α12[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]) -
            sum(sm.bigM*mas.ujk[j,k,m]*subp.α13[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]) -
            sum(sm.bigM*mas.ukl[k,l,m]*subp.α14[k,l,m] for k=1:sm.K for l=1:sm.L for m=1:sm.Mkl[k,l])
            )
        MOI.submit(mas.m, MOI.LazyConstraint(cb_data), cut); push!(tmpnf,1)
        # push!(fcuts, (α5=subp.α5,α6=subp.α6,α7=subp.α7,α8=subp.α8,α9=subp.α9,α10=subp.α10,α12=subp.α12,α13=subp.α13,α14=subp.α14))            

    end
    push!(tmpit,1)
    return 
end

function ClassicBD(sm,wdg,rort,buf)#path,ins,
    # lex1 = build_scndmip(1,0,sm,1e-4);
    # optimize!(lex1)
    # f1min,f2max = getobjval(lex1,sm)
    # lex2 = build_scndmip(0,1,sm,1e-4);
    # optimize!(lex2)
    # @show f1max,f2min = getobjval(lex2,sm)

    mip = smSCND_MIP(sm,1e-4);
    vSolve(mip, Inf,0,method=:lex); #vSolve(mip, Inf,method=:lex);
    f1,f2 = getY_N(mip)

    ##Benders Decomposition Algorithm 
    # solpt,tb2= CBD(sm,wdg,rort,f1min,f2max,f1max,f2min);
    solpt,tb,tb2= CBD(sm,wdg,rort,buf,round(f1[1]),round(f1[2]),round(f2[1]),round(f2[2])); #f1[1],f1[2],f2[1],f2[2]
    # solpt,tb,tb2= CBD_fixedW(sm,wdg,rort,buf,round(f1[1]),round(f1[2]),round(f2[1]),round(f2[2]),[0.85798,0.93959,0.805,0.75515, 0.83439]); #f1[1],f1[2],f2[1],f2[2]
    # ndp = NDfilter(solpt)
    return solpt,tb,tb2
end
function ClassicBD_trunc(sm,wdg,rort,buf)#path,ins,
    mip = smSCND_MIP(sm,1e-4);
    vSolve(mip, Inf,0,method=:lex); #vSolve(mip, Inf,method=:lex);
    f1,f2 = getY_N(mip)
    # solpt,tb,tb2= CBD_trunc(sm,wdg,rort,buf,round(f1[1]),round(f1[2]),round(f2[1]),round(f2[2])); #f1[1],f1[2],f2[1],f2[2]
    solpt,tb,tb2= CBD_trunc_fixedW(sm,wdg,rort,buf,round(f1[1]),round(f1[2]),round(f2[1]),round(f2[2]),at.w); #[0.8579,0.9395,0.805,0.7551, 0.8343][0.85798,0.93959,0.805,0.75515, 0.83439]
    return solpt,tb,tb2
end
function CBD(sm,wdg,rort,buf,f1min,f2max,f1max,f2min)
    # global fcuts = []; global ocuts = [] ;
    tb = DataFrame(fname = [], w = [], obj1 = [], obj2 = [], nfcut = [], nocut = [], subCPU = [] ,MPCPU=[], Added = []);
    tb2 = DataFrame( pt1 = [], pt2 = [], w = [], lb = [], val = [], obj1 = [] ,obj2=[], wobj=[], Added = []);
    L = [];  global subtime,tmpnf,tmpno,tmpit = [.0],[0],[0],[0]; global YN = [[f1min,f2max],[f1max,f2min]];
    λ1 = abs(f2max-f2min);  λ2 = abs(f1max-f1min); ratio = λ2/λ1; 
    if rort == "r"
        global w = round(1/(1+ratio); digits=wdg)
    else
        global w = trunc(1/(1+ratio); digits=wdg)
    end
    # println("w = ", w)
    global mas = Master(w,sm,1e-4); global dsp = DualSubP(w,1-w,sm,1e-4); MOI.set(mas.m, MOI.LazyConstraintCallback(), Classic_lazy_callback);
    tCPU = @CPUelapsed optimize!(mas.m); yt_1,yt_2 = getobjdual(mas.m,dsp.m,sm); 
    lb = w*f1min + (1-w)*f2max;  val = w*yt_1 + (1-w)*yt_2
    # if ( val + 10 < lb) #&& f1min < yt_1 < f1max && f2min < yt_2 < f2max
    push!(YN, [yt_1,yt_2]); push!(L,[[yt_1,yt_2],(f1max,f2min)],[(f1min,f2max),[yt_1,yt_2]])
    # println(yt_1," ",yt_2)
    # end
    push!(tb, (sm.insfile[end-11:end-4], w, yt_1, yt_2, sum(tmpnf), sum(tmpno), sum(subtime), tCPU, "Y"))
    push!(tb2, [f1min,f1max,w,lb,val, yt_1, yt_2,objective_value(mas.m),"Y"]) #insert!(tb, 1, ("$fname",w, objective_value(mas.m), yt_1, yt_2, sum(tmpnf), sum(tmpno), sum(tmpit), sum(subtime), tCPU) )
    while L != []
        insert = "N"
        global (subtime,tmpnf,tmpno,tmpit) = [.0],[0],[0],[0]; 
        pt = pop!(L)
        λ1 = abs(pt[1][2]-pt[2][2])
        λ2 = abs(pt[2][1]-pt[1][1])
        ratio = λ2/λ1
        if rort == "r"
            w = round(1/(1+ratio); digits=wdg)
        else
            w = trunc(1/(1+ratio); digits=wdg)
        end
        # println("w = ", w)
        if λ1 ==0 || w ∈ tb.w || round(w;digits=3)>=.99 || trunc(w;digits=2) <= .01 || w == NaN  
        # if λ1 ==0 || w ∈ tb.w || round(w;digits=3)>=.95 || trunc(w;digits=2) <= .01 || w == NaN  
            @goto Nextiter 
        end
        lb = w*pt[1][1] + (1-w)*pt[1][2]

        global mas = Master(w,sm,1e-4); global dsp = DualSubP(w,1-w,sm,1e-4); 
        MOI.set(mas.m, MOI.LazyConstraintCallback(), Classic_lazy_callback);
        tCPU = @CPUelapsed optimize!(mas.m)
        if termination_status(mas.m) in [MOI.OPTIMAL,MOI.LOCALLY_SOLVED]
            yt_1,yt_2 = getobjdual(mas.m,dsp.m,sm) 
            val = objective_value(mas.m) #w*yt_1 + (1-w)*yt_2
            if (val + buf < lb) && pt[1][1]< yt_1 < pt[2][1] && pt[2][2] < yt_2 < pt[1][2]#dominated([yt_1,yt_2],YN)==false &&
                println(yt_1," ",yt_2)
                insert = "Y"
                push!(YN, [yt_1,yt_2]); push!(L,[pt[1],[yt_1,yt_2]],[[yt_1,yt_2],pt[2]])
            end    
            push!(tb, (sm.insfile[end-11:end-4], w, yt_1, yt_2, sum(tmpnf), sum(tmpno), sum(subtime), tCPU, insert)) # push!(tb, (sm.insfile[end-11:end-4],, w, lb,val, objective_value(mas.m),yt_1,yt_2, sum(tmpnf), sum(tmpno), sum(tmpit),sum(subtime), tCPU))
            push!(tb2, [pt[1][1],pt[2][1],w,lb,val,yt_1,yt_2,objective_value(mas.m),insert])
            # insert!(tb, size(tb,1), ("$fname",w, objective_value(mas.m),yt_1,yt_2, sum(tmpnf), sum(tmpno), sum(tmpit),round(sum(subtime); digits=2), round(tCPU; digits=2)))
        end    
        @label Nextiter
    end

    return YN,tb,tb2 
end
function CBD_trunc(sm,wdg,rort,buf,f1min,f2max,f1max,f2min)
    # global fcuts = []; global ocuts = [] ;
    tb = DataFrame(fname = [], w = [], obj1 = [], obj2 = [], nfcut = [], nocut = [], subCPU = [] ,MPCPU=[], Added = []);
    tb2 = DataFrame( pt1 = [], pt2 = [], w = [], lb = [], val = [], obj1 = [] ,obj2=[], mobj = [], Added = []);
    L = [];  global subtime,tmpnf,tmpno,tmpit = [.0],[0],[0],[0]; global YN = [[f1min,f2max],[f1max,f2min]];
    λ1 = abs(f2max-f2min);  λ2 = abs(f1max-f1min); ratio = λ2/λ1; 
    if rort == "r"
        global w = round(1/(1+ratio); digits=wdg)
    else
        global w = trunc(1/(1+ratio); digits=wdg)
    end    
    # println("w = ", w)
    global mas = Master(w,sm,1e-4); global dsp = DualSubP(w,1-w,sm,1e-4); MOI.set(mas.m, MOI.LazyConstraintCallback(), Classic_lazy_callback_trunc);
    tCPU = @CPUelapsed optimize!(mas.m); yt_1,yt_2 = getobjdual(mas.m,dsp.m,sm); 
    lb = w*f1min + (1-w)*f2max;  val = w*yt_1 + (1-w)*yt_2
    push!(YN, [yt_1,yt_2]); push!(L,[[yt_1,yt_2],(f1max,f2min)],[(f1min,f2max),[yt_1,yt_2]])
    push!(tb, (sm.insfile[end-11:end-4], w, yt_1, yt_2, sum(tmpnf), sum(tmpno), sum(subtime), tCPU, "Y"))
    push!(tb2, [f1min,f1max,w,lb,val, yt_1, yt_2, objective_value(mas.m),"Y"]) #insert!(tb, 1, ("$fname",w, objective_value(mas.m), yt_1, yt_2, sum(tmpnf), sum(tmpno), sum(tmpit), sum(subtime), tCPU) )
    while L != []
        insert = "N"
        global (subtime,tmpnf,tmpno,tmpit) = [.0],[0],[0],[0]; 
        pt = pop!(L)
        λ1 = abs(pt[1][2]-pt[2][2])
        λ2 = abs(pt[2][1]-pt[1][1])
        ratio = λ2/λ1
        if rort == "r"
            w = round(1/(1+ratio); digits=wdg)
        else
            w = trunc(1/(1+ratio); digits=wdg)
        end
        # println("w = ", w)
        if λ1 ==0 || w ∈ tb.w || round(w;digits=3)>=.991 || trunc(w;digits=2) <= .01 || w == NaN  
            @goto Nextiter 
        end
        lb = w*pt[1][1] + (1-w)*pt[1][2]

        global mas = Master(w,sm,1e-4); global dsp = DualSubP(w,1-w,sm,1e-4); 
        MOI.set(mas.m, MOI.LazyConstraintCallback(), Classic_lazy_callback_trunc);
        tCPU = @CPUelapsed optimize!(mas.m)
        if termination_status(mas.m) in [MOI.OPTIMAL,MOI.LOCALLY_SOLVED]
            yt_1,yt_2 = getobjdual(mas.m,dsp.m,sm) 
            val = objective_value(mas.m) #w*yt_1 + (1-w)*yt_2
            
            if (val + buf < lb) && pt[1][1]< yt_1 < pt[2][1] && pt[2][2] < yt_2 < pt[1][2]#dominated([yt_1,yt_2],YN)==false &&
                println(yt_1," ",yt_2)
                insert = "Y"
                push!(YN, [yt_1,yt_2]); push!(L,[pt[1],[yt_1,yt_2]],[[yt_1,yt_2],pt[2]])
            end    
            push!(tb, (sm.insfile[end-11:end-4], w, yt_1, yt_2, sum(tmpnf), sum(tmpno), sum(subtime), tCPU, insert)) # push!(tb, (sm.insfile[end-11:end-4],, w, lb,val, objective_value(mas.m),yt_1,yt_2, sum(tmpnf), sum(tmpno), sum(tmpit),sum(subtime), tCPU))
            push!(tb2, [pt[1][1],pt[2][1],w,lb,w*yt_1 + (1-w)*yt_2,yt_1,yt_2,objective_value(mas.m),insert])
        end    
        @label Nextiter
    end

    return YN,tb,tb2 
end

#################Accelerated Benders
function Acc_lazy_callback(cb_data)
    yb = callback_value.(cb_data, mas.y);    ubij = callback_value.(cb_data, mas.uij);    ubjk = callback_value.(cb_data, mas.ujk)
    ubkl = callback_value.(cb_data, mas.ukl);    θb = callback_value(cb_data, mas.θ)
    # global cbint+=1
    subp = solve_subp(dsp,yb,ubij,ubjk,ubkl);    
    if subp.res == :OptimalityCut
        
        # @info "Optimality cut found"
        if θb ≥ round(subp.obj; digits=num)
        # if θb ≥ trunc(subp.obj; digits=num)
            return
        else
            head0 = Vector{Cint}(undef, numcon); bx = Vector{Cdouble}(undef, numcon)
            CPLEX.CPXgetbhead(backend(dsp.m).env, backend(dsp.m).lp, head0, bx); head = head0.+1;
            c1 = copy(rhs1); c2 = copy(rhs2); 
            CPLEX.CPXftran(backend(dsp.m).env, backend(dsp.m).lp, c1)
            CPLEX.CPXftran(backend(dsp.m).env, backend(dsp.m).lp, c2)
            pi1,pi2 = value.(all_variables(dsp.m)), value.(all_variables(dsp.m))
            for (idx,i) in enumerate(head)
                if i > 0
                    pi1[i] = c1[idx]
                    pi2[i] = c2[idx]
                end
            end
            ppi1,ppi2 = copy(pi1),copy(pi2) 
            if [ppi1,ppi2] ∉ dfo.ppi
                push!(dfo, [w, [ppi1,ppi2] ,subp.α5,subp.α6])
            end
            # push!(dfo, [w,round.(copy(pi1); digits=1),round.(copy(pi2); digits=1),subp.α5,subp.α6])
            # push!(dfo, [w,copy(pi1)[2*(sm.J+sm.K)+sm.L+sm.I+1:end],copy(pi2)[2*(sm.J+sm.K)+sm.L+sm.I+1:end],subp.α5,subp.α6]);
            # push!(ocuts, (obj=subp.obj,α5=subp.α5,α6=subp.α6,α7=subp.α7,α8=subp.α8,α9=subp.α9,α10=subp.α10,α12=subp.α12,α13=subp.α13,α14=subp.α14))
            cut = @build_constraint( mas.θ ≥ 
                dot(sm.demand,subp.α5)-dot(sm.cas,subp.α6)-
                sum(sm.cap[j]*subp.α7[j,t]*mas.y[j,t] for j=1:sm.J for t=1:2) -
                sum(sm.cad[k]*subp.α8[k,t]*mas.y[k+sm.J,t] for k=1:sm.K for t=1:2) +
                sum(sm.Lcapasp[i][j][m]*subp.α9[i,j,m]*mas.uij[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]) +
                sum(sm.Lcapapd[j][k][m]*subp.α10[j,k,m]*mas.ujk[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]) -
                sum(sm.bigM*subp.α12[i,j,m]*mas.uij[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]) -
                sum(sm.bigM*subp.α13[j,k,m]*mas.ujk[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]) -
                sum(sm.bigM*subp.α14[k,l,m]*mas.ukl[k,l,m] for k=1:sm.K for l=1:sm.L for m=1:sm.Mkl[k,l])
            )
            MOI.submit(mas.m, MOI.LazyConstraint(cb_data), cut)
            push!(tmpno,1)
        end
    # else #
    elseif subp.res == :FeasibilityCut
        # @info "Feasibility cut found"
        cut = @build_constraint( 0 ≥ dot(sm.demand,subp.α5)-dot(sm.cas,subp.α6)-
            sum(sm.cap[j]*subp.α7[j,t]*mas.y[j,t] for j=1:sm.J for t=1:2) -
            sum(sm.cad[k]*subp.α8[k,t]*mas.y[k+sm.J,t] for k=1:sm.K for t=1:2) +
            sum(sm.Lcapasp[i][j][m]*subp.α9[i,j,m]*mas.uij[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]) +
            sum(sm.Lcapapd[j][k][m]*subp.α10[j,k,m]*mas.ujk[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]) -
            sum(sm.bigM*mas.uij[i,j,m]*subp.α12[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]) -
            sum(sm.bigM*mas.ujk[j,k,m]*subp.α13[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]) -
            sum(sm.bigM*mas.ukl[k,l,m]*subp.α14[k,l,m] for k=1:sm.K for l=1:sm.L for m=1:sm.Mkl[k,l]))
        MOI.submit(mas.m, MOI.LazyConstraint(cb_data), cut)
        # push!(dff, [w, (α5=subp.α5,α6=subp.α6,α7=subp.α7,α8=subp.α8,α9=subp.α9,α10=subp.α10,α12=subp.α12,α13=subp.α13,α14=subp.α14)])
        # push!(fcuts, (α5=subp.α5,α6=subp.α6,α7=subp.α7,α8=subp.α8,α9=subp.α9,α10=subp.α10,α12=subp.α12,α13=subp.α13,α14=subp.α14))
        push!(tmpnf,1)
    end
    push!(tmpit,1)
    return 
end
function Acc_lazy_callback_trunc(cb_data)
    yb = callback_value.(cb_data, mas.y);    ubij = callback_value.(cb_data, mas.uij);    ubjk = callback_value.(cb_data, mas.ujk)
    ubkl = callback_value.(cb_data, mas.ukl);    θb = callback_value(cb_data, mas.θ)
    # global cbint+=1
    subp = solve_subp(dsp,yb,ubij,ubjk,ubkl);    
    if subp.res == :OptimalityCut
        if θb ≥ trunc(subp.obj; digits=num)
            return
        else
            head0 = Vector{Cint}(undef, numcon); bx = Vector{Cdouble}(undef, numcon)
            CPLEX.CPXgetbhead(backend(dsp.m).env, backend(dsp.m).lp, head0, bx); head = head0.+1;
            c1 = copy(rhs1); c2 = copy(rhs2); 
            CPLEX.CPXftran(backend(dsp.m).env, backend(dsp.m).lp, c1)
            CPLEX.CPXftran(backend(dsp.m).env, backend(dsp.m).lp, c2)
            pi1,pi2 = value.(all_variables(dsp.m)), value.(all_variables(dsp.m))
            for (idx,i) in enumerate(head)
                if i > 0
                    pi1[i] = c1[idx]
                    pi2[i] = c2[idx]
                end
            end
           
            ppi1,ppi2 = copy(pi1),copy(pi2) 
            if [ppi1,ppi2] ∉ dfo.ppi
                push!(dfo, [w, [ppi1,ppi2] ,subp.α5,subp.α6])
            end
           
            cut = @build_constraint( mas.θ ≥ 
                dot(sm.demand,subp.α5)-dot(sm.cas,subp.α6)-
                sum(sm.cap[j]*subp.α7[j,t]*mas.y[j,t] for j=1:sm.J for t=1:2) -
                sum(sm.cad[k]*subp.α8[k,t]*mas.y[k+sm.J,t] for k=1:sm.K for t=1:2) +
                sum(sm.Lcapasp[i][j][m]*subp.α9[i,j,m]*mas.uij[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]) +
                sum(sm.Lcapapd[j][k][m]*subp.α10[j,k,m]*mas.ujk[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]) -
                sum(sm.bigM*subp.α12[i,j,m]*mas.uij[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]) -
                sum(sm.bigM*subp.α13[j,k,m]*mas.ujk[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]) -
                sum(sm.bigM*subp.α14[k,l,m]*mas.ukl[k,l,m] for k=1:sm.K for l=1:sm.L for m=1:sm.Mkl[k,l])
            )
            MOI.submit(mas.m, MOI.LazyConstraint(cb_data), cut)
            push!(tmpno,1)
        end
    # else #
    elseif subp.res == :FeasibilityCut
        # @info "Feasibility cut found"
        cut = @build_constraint( 0 ≥ dot(sm.demand,subp.α5)-dot(sm.cas,subp.α6)-
            sum(sm.cap[j]*subp.α7[j,t]*mas.y[j,t] for j=1:sm.J for t=1:2) -
            sum(sm.cad[k]*subp.α8[k,t]*mas.y[k+sm.J,t] for k=1:sm.K for t=1:2) +
            sum(sm.Lcapasp[i][j][m]*subp.α9[i,j,m]*mas.uij[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]) +
            sum(sm.Lcapapd[j][k][m]*subp.α10[j,k,m]*mas.ujk[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]) -
            sum(sm.bigM*mas.uij[i,j,m]*subp.α12[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]) -
            sum(sm.bigM*mas.ujk[j,k,m]*subp.α13[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]) -
            sum(sm.bigM*mas.ukl[k,l,m]*subp.α14[k,l,m] for k=1:sm.K for l=1:sm.L for m=1:sm.Mkl[k,l]))
        MOI.submit(mas.m, MOI.LazyConstraint(cb_data), cut)
        # push!(dff, [w, (α5=subp.α5,α6=subp.α6,α7=subp.α7,α8=subp.α8,α9=subp.α9,α10=subp.α10,α12=subp.α12,α13=subp.α13,α14=subp.α14)])
        # push!(fcuts, (α5=subp.α5,α6=subp.α6,α7=subp.α7,α8=subp.α8,α9=subp.α9,α10=subp.α10,α12=subp.α12,α13=subp.α13,α14=subp.α14))
        push!(tmpnf,1)
    end
    push!(tmpit,1)
    return 
end
function AcceleratedBD(sm,pcut,wdg,rort,dg,buf)
    # m-dimensional model
    global rhs1 = [[(sm.vcs[i]+sm.tcp[i][j][m]) for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]];[sm.tcd[j][k][m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]];[sm.tcc[k][l][m] for k=1:sm.K for l=1:sm.L for m=1:sm.Mkl[k,l]]; [sm.vcp[j,t] for t=1:2 for j=1:sm.J ];[sm.vcd[k,t] for t=1:2 for k=1:sm.K ]]; 
    global rhs2 = [[(sm.ves[i]+sm.cep[i][j][m]) for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]];[sm.ced[j][k][m] for j=1:sm.J for k=1:sm.K for m=1:1 for m=1:sm.Mjk[j,k]];[sm.cec[k][l][m] for k=1:sm.K for l=1:sm.L for m=1:sm.Mkl[k,l]]; [sm.vep[j,t] for t=1:2 for j=1:sm.J ];[sm.ved[k,t] for t=1:2 for k=1:sm.K ]];
    global len = [length(dsp.α1),length(dsp.α2),length(dsp.α3),length(dsp.α4),length(dsp.α5),length(dsp.α6),length(dsp.α7),length(dsp.α8),length(dsp.α9),length(dsp.α10),length(dsp.α12),length(dsp.α13),length(dsp.α14)]
    global numcon = num_constraints(dsp.m, AffExpr, MOI.LessThan{Float64})
    global objc = [-[sm.cap[j] for t=1:2 for j=1:sm.J];-[sm.cad[k] for t=1:2 for k=1:sm.K];[sm.Lcapasp[i][j][m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]];[sm.Lcapapd[j][k][m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]];-sm.bigM*ones(sum(sm.Mij)+sum(sm.Mjk)+sum(sm.Mkl))]
    
    mip = smSCND_MIP(sm,1e-4);
    vSolve(mip,Inf,0,method=:lex );#vSolve(mip,Inf,method=:lex );
    f1,f2 = getY_N(mip)
    ###Accelerated Benders Decomposition Algorithm 
    solpt,tb,tb2 = AccBD(sm,pcut,wdg,rort,dg,buf,round(f1[1]),round(f1[2]),round(f2[1]),round(f2[2]))
    # solpt,tb,tb2 = AccBD_fixedW(sm,pcut,wdg,rort,dg,buf,round(f1[1]),round(f1[2]),round(f2[1]),round(f2[2]),[0.85798,0.93959,0.805,0.75515, 0.83439])
    return solpt,tb,tb2
end
function AcceleratedBD_trunc(sm,pcut,wdg,rort,dg,buf)
    global rhs1 = [[(sm.vcs[i]+sm.tcp[i][j][m]) for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]];[sm.tcd[j][k][m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]];[sm.tcc[k][l][m] for k=1:sm.K for l=1:sm.L for m=1:sm.Mkl[k,l]]; [sm.vcp[j,t] for t=1:2 for j=1:sm.J ];[sm.vcd[k,t] for t=1:2 for k=1:sm.K ]]; 
    global rhs2 = [[(sm.ves[i]+sm.cep[i][j][m]) for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]];[sm.ced[j][k][m] for j=1:sm.J for k=1:sm.K for m=1:1 for m=1:sm.Mjk[j,k]];[sm.cec[k][l][m] for k=1:sm.K for l=1:sm.L for m=1:sm.Mkl[k,l]]; [sm.vep[j,t] for t=1:2 for j=1:sm.J ];[sm.ved[k,t] for t=1:2 for k=1:sm.K ]];
    global len = [length(dsp.α1),length(dsp.α2),length(dsp.α3),length(dsp.α4),length(dsp.α5),length(dsp.α6),length(dsp.α7),length(dsp.α8),length(dsp.α9),length(dsp.α10),length(dsp.α12),length(dsp.α13),length(dsp.α14)]
    global numcon = num_constraints(dsp.m, AffExpr, MOI.LessThan{Float64})
    global objc = [-[sm.cap[j] for t=1:2 for j=1:sm.J];-[sm.cad[k] for t=1:2 for k=1:sm.K];[sm.Lcapasp[i][j][m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]];[sm.Lcapapd[j][k][m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]];-sm.bigM*ones(sum(sm.Mij)+sum(sm.Mjk)+sum(sm.Mkl))]
    mip = smSCND_MIP(sm,1e-4);
    vSolve(mip,Inf,0,method=:lex );#vSolve(mip,Inf,method=:lex );
    f1,f2 = getY_N(mip)
    # solpt,tb,tb2 = AccBD_trunc(sm,pcut,wdg,rort,dg,buf,round(f1[1]),round(f1[2]),round(f2[1]),round(f2[2]))
    solpt,tb,tb2 = AccBD_trunc_fixedW(sm,pcut,wdg,rort,dg,buf,round(f1[1]),round(f1[2]),round(f2[1]),round(f2[2]),at.w)#[0.8579,0.9395,0.805,0.7551, 0.8343]
    return solpt,tb,tb2#,cutw
end
function AccBD(sm,pcut,wdg,rort,dg,buf,f1min,f2max,f1max,f2min)
    tb = DataFrame(fname = [], w = [], obj1 = [], obj2 = [], preOcut = [], nfcut = [], nocut = [], cutCPU =[], subCPU = [] ,MPCPU=[], Added = []); 
    tb2 = DataFrame( pt1 = [], pt2 = [], w = [], lb = [], val = [], obj1 = [] ,obj2=[], mobj=[], Added = []);
    global dfo = DataFrame(w=[], ppi=[], α5=[],α6=[]); global dff = DataFrame(w=[],cuts=[])
    L = []; global (subtime,tmpnf,tmpno) = ([.0],[0],[0]); global YN = [[f1min,f2max],[f1max,f2min]]; 
    λ1 = abs(f2max-f2min);  λ2 = abs(f1max-f1min); ratio = λ2/λ1; 
    if rort == "r"
        global w = round(1/(1+ratio); digits=wdg)
    else
        global w = trunc(1/(1+ratio); digits=wdg)
    end
    # println("w = ", w)
    global mas = Master(w,sm,1e-4); global dsp = DualSubP(w,1-w,sm,1e-4); MOI.set(mas.m, MOI.LazyConstraintCallback(), Acc_lazy_callback);
    tCPU = @CPUelapsed optimize!(mas.m); yt_1,yt_2 = getobjdual(mas.m,dsp.m,sm); 
    lb = w*f1min + (1-w)*f2max;  val = w*yt_1 + (1-w)*yt_2

    push!(YN, [yt_1,yt_2])# println(yt_1," ",yt_2)
    push!(L,[[yt_1,yt_2],(f1max,f2min)],[(f1min,f2max),[yt_1,yt_2]])
    push!(tb, (sm.insfile[end-11:end-4], w, yt_1, yt_2, 0, sum(tmpnf), sum(tmpno), 0, sum(subtime), tCPU, "Y"))
    push!(tb2, [f1min,f1max,w,lb,val,yt_1,yt_2,objective_value(mas.m),"Y"])
    # cutw = DataFrame(ins = [], current=[], from=[], idx=[])
    while L != []
        insert = "N" 
        subtime,tmpnf,tmpno = [.0],[0],[0]; 
        pt = pop!(L)
        λ1 = abs(pt[1][2]-pt[2][2])
        λ2 = abs(pt[2][1]-pt[1][1])
        ratio = λ2/λ1
        if rort == "r"
            w = round(1/(1+ratio); digits=wdg)
        else
            w = trunc(1/(1+ratio); digits=wdg)
        end
        # println(" w = ", w) 
        if λ1 ==0 || w ∈ tb.w || round(w;digits=3)>=.991 || trunc(w;digits=2) <= .01 || w == NaN  
        # if λ1 ==0 || w ∈ tb.w || round(w;digits=3)>=.98 || trunc(w;digits=2) <= .01 || w == NaN  
            @goto Nextiter 
        end
        lb = w*pt[1][1] + (1-w)*pt[1][2]

        # Add only violated cuts  
        oid = []; ogap = DataFrame(id=[], gap=[]); #fid = []
        mpobj = value(mas.θ); 
        vy = value.(mas.m[:y]); vuij = value.(mas.m[:uij]); vujk = value.(mas.m[:ujk]); vukl = value.(mas.m[:ukl]);
        masv = [[value.(mas.m[:y])[j,t] for t=1:2 for j=1:sm.J];[value.(mas.m[:y])[k,t] for t=1:2 for k=sm.J+1:sm.J+sm.K];
            [value.(mas.m[:uij])[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]];
            [value.(mas.m[:ujk])[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]];
            [value.(mas.m[:uij])[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]];
            [value.(mas.m[:ujk])[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]];
            [value.(mas.m[:ukl])[k,l,m] for k=1:sm.K for l=1:sm.L for m=1:sm.Mkl[k,l]]]

        t0 = time()       
        global mas = Master(w,sm,1e-4); global dsp = DualSubP(w,1-w,sm,1e-4); MOI.set(mas.m, MOI.LazyConstraintCallback(), Acc_lazy_callback);
        mvar = [[mas.y[j,t] for j=1:sm.J+sm.K for t=1:2]; [mas.uij[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]]; [mas.ujk[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]]; [mas.uij[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]]; [mas.ujk[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]];[mas.ukl[k,l,m] for k=1:sm.K for l=1:sm.L for m=1:sm.Mkl[k,l]]]

        for o=1:size(dfo,1)
            ov = dot(sm.demand,dfo.α5[o])-dot(sm.cas,dfo.α6[o])+dot((objc.*trunc.(w*dfo.ppi[o][1][1+sum(len[1:6]):end]+(1-w)*dfo.ppi[o][2][1+sum(len[1:6]):end]; digits=dg)), masv)
            if ov > mpobj && all(i->i>=0,(w*dfo.ppi[o][1]+(1-w)*dfo.ppi[o][2])) == true            
                # @constraint(mas.m, dot(sm.demand,dfo.α5[o])-dot(sm.cas,dfo.α6[o])+dot((objc.*trunc.(w*dfo.ppi[o][1][1+sum(len[1:6]):end]+(1-w)*dfo.ppi[o][2][1+sum(len[1:6]):end]; digits=dg)), mvar)  <= mas.θ)                
                # println("optimality cut = ",ov," current θ =", value(mas.θ))
                push!(oid,o); push!(ogap, (o, abs(ov-mpobj)) )
            # else
                # @show all(i->i>=0,(w*dfo.ppi[o][1]+(1-w)*dfo.ppi[o][2]))
            end
        end
        # println(oid)
        
        ## random selection
        # pick = min(pcut,length(oid))
        # ocutid = sample(oid, pick; replace=false)
        # for o in ocutid
        #     @constraint(mas.m, dot(sm.demand,dfo.α5[o])-dot(sm.cas,dfo.α6[o])+dot((objc.*trunc.(w*dfo.ppi[o][1][1+sum(len[1:6]):end]+(1-w)*dfo.ppi[o][2][1+sum(len[1:6]):end]; digits=dg)), mvar)  <= mas.θ)
        # end

        # select from sorted violated cuts
        pick = min(pcut,length(oid))
        # sort!(ogap, [:gap], rev=true) #most violated
        sort!(ogap, [:gap]) #least violated
        ocutid = ogap.id[1:pick]
        for o in ocutid
            # push!(cutw, [sm.insfile[end-11:end-4], w, dfo.w[o], o])
            @constraint(mas.m, dot((objc.*trunc.(w*dfo.ppi[o][1][1+sum(len[1:6]):end]+(1-w)*dfo.ppi[o][2][1+sum(len[1:6]):end]; digits=dg)), mvar)  <= mas.θ)
        end
        
        cutime = round(time()-t0; digits=2)
        tCPU = @CPUelapsed optimize!(mas.m)
        if termination_status(mas.m) in [MOI.OPTIMAL,MOI.LOCALLY_SOLVED]
            yt_1,yt_2 = getobjdual(mas.m,dsp.m,sm) 
            # @show termination_status(mas.m)
            val = w*yt_1 + (1-w)*yt_2
            
            if (val + buf < lb) && pt[1][1]< yt_1 < pt[2][1] && pt[2][2] < yt_2 < pt[1][2]#round.([yt_1,yt_2]) ∉ YN && val + 10 < lb #&& f1min < yt_1 && yt_2 < f2max#[yt_1,yt_2] ∉ YN #
                println(yt_1," ",yt_2)
                insert = "Y"
                push!(YN, [yt_1,yt_2]); push!(L,[pt[1],[yt_1,yt_2]],[[yt_1,yt_2],pt[2]])
                # push!(YN,round.([yt_1,yt_2]; digits=2)); push!(L,[pt[1],round.([yt_1,yt_2]; digits=2)],[round.([yt_1,yt_2]; digits=2),pt[2]])
            end 
            #ocutid = oid
            push!(tb, (sm.insfile[end-11:end-4], w, yt_1, yt_2, length(ocutid),sum(tmpnf), sum(tmpno), cutime, sum(subtime), tCPU, insert)) 
            push!(tb2, [pt[1][1],pt[2][1],w,lb,val,yt_1,yt_2,objective_value(mas.m),insert]) 
        end   
        @label Nextiter
    end
    return YN,tb,tb2##,cutw
end
function AccBD_trunc(sm,pcut,wdg,dgrort,dg,buf,f1min,f2max,f1max,f2min)
    tb = DataFrame(fname = [], w = [], obj1 = [], obj2 = [], preOcut = [], nfcut = [], nocut = [], cutCPU =[], subCPU = [] ,MPCPU=[], Added = []);   
    tb2 = DataFrame( pt1 = [], pt2 = [], w = [], lb = [], val = [], obj1 = [] ,obj2=[], mobj=[], Added = []);
    global dfo = DataFrame(w=[], ppi=[], α5=[],α6=[]); global dff = DataFrame(w=[],cuts=[])# global (pi1,pi2,fcuts,ocuts) = [],[],[],[] 
    L = []; global (subtime,tmpnf,tmpno,tmpit) = ([.0],[0],[0],[0]); global YN = [[f1min,f2max],[f1max,f2min]]; # global ocuts = []; 
    λ1 = abs(f2max-f2min);  λ2 = abs(f1max-f1min); ratio = λ2/λ1; 
    if rort == "r"
        global w = round(1/(1+ratio); digits=wdg)
    else
        global w = trunc(1/(1+ratio); digits=wdg)
    end
    # println("w = ", w);  
    global mas = Master(w,sm,1e-4); global dsp = DualSubP(w,1-w,sm,1e-4); MOI.set(mas.m, MOI.LazyConstraintCallback(), Acc_lazy_callback_trunc);
    tCPU = @CPUelapsed optimize!(mas.m); yt_1,yt_2 = getobjdual(mas.m,dsp.m,sm); 
    lb = w*f1min + (1-w)*f2max;  val = w*yt_1 + (1-w)*yt_2

    push!(YN, [yt_1,yt_2])# println(yt_1," ",yt_2)
    push!(L,[[yt_1,yt_2],(f1max,f2min)],[(f1min,f2max),[yt_1,yt_2]])
    push!(tb, (sm.insfile[end-11:end-4], w, yt_1, yt_2, 0, sum(tmpnf), sum(tmpno), 0, sum(subtime), tCPU, "Y"))
    push!(tb2, [f1min,f1max,w,lb,val,yt_1,yt_2,objective_value(mas.m),"Y"])
    cutw = DataFrame(ins = [], current=[], from=[], idx=[])
    while L != []
        insert = "N" 
        subtime,tmpnf,tmpno,tmpit = [.0],[0],[0],[0]; 
        pt = pop!(L)
        λ1 = abs(pt[1][2]-pt[2][2])
        λ2 = abs(pt[2][1]-pt[1][1])
        ratio = λ2/λ1
        if rort == "r"
            w = round(1/(1+ratio); digits=wdg)
        else
            w = trunc(1/(1+ratio); digits=wdg)
        end
        # println(" w = ", w) 
        if λ1 ==0 || w ∈ tb.w || round(w;digits=3)>=.99 || trunc(w;digits=2) <= .01 || w == NaN
            @goto Nextiter
        end
        lb = w*pt[1][1] + (1-w)*pt[1][2]

        # Add only violated cuts
        oid = []; ogap = DataFrame(id=[], gap=[])    #fid = []
        mpobj = value(mas.θ);
        masv = [[value.(mas.m[:y])[j,t] for t=1:2 for j=1:sm.J];[value.(mas.m[:y])[k,t] for t=1:2 for k=sm.J+1:sm.J+sm.K];
            [value.(mas.m[:uij])[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]];
            [value.(mas.m[:ujk])[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]];
            [value.(mas.m[:uij])[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]];
            [value.(mas.m[:ujk])[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]];
            [value.(mas.m[:ukl])[k,l,m] for k=1:sm.K for l=1:sm.L for m=1:sm.Mkl[k,l]]]

        t0 = time()
        
        for o=1:size(dfo,1)
            ov = dot(sm.demand,dfo.α5[o])-dot(sm.cas,dfo.α6[o])+dot((objc.*trunc.(w*dfo.ppi[o][1][1+sum(len[1:6]):end]+(1-w)*dfo.ppi[o][2][1+sum(len[1:6]):end]; digits=dg)), masv)
            # ov = dot(sm.demand,dfo.α5[o])-dot(sm.cas,dfo.α6[o])+dot((objc.*round.(w*dfo.ppi[o][1][1+sum(len[1:6]):end]+(1-w)*dfo.ppi[o][2][1+sum(len[1:6]):end]; digits=dg)), masv)
            if ov > mpobj && all(i->i>=0,(w*dfo.ppi[o][1]+(1-w)*dfo.ppi[o][2])) == true
                # @constraint(mas.m, dot(sm.demand,dfo.α5[o])-dot(sm.cas,dfo.α6[o])+dot((objc.*trunc.(w*dfo.ppi[o][1][1+sum(len[1:6]):end]+(1-w)*dfo.ppi[o][2][1+sum(len[1:6]):end]; digits=dg)), mvar)  <= mas.θ)
                push!(oid,o); push!(ogap, (o, abs(ov-mpobj)) )
            end
        end

        # random selection
        # pick = min(pcut,length(oid))
        # ocutid = sample(oid, pick; replace=false)
        # for o in ocutid
        #     @constraint(mas.m, dot(sm.demand,dfo.α5[o])-dot(sm.cas,dfo.α6[o])+dot((objc.*trunc.(w*dfo.ppi[o][1][1+sum(len[1:6]):end]+(1-w)*dfo.ppi[o][2][1+sum(len[1:6]):end]; digits=dg)), mvar)  <= mas.θ)
        # end

        #select from sorted violated cuts
        pick = min(pcut,length(oid))
        # sort!(ogap, [:gap], rev=true) #most violated
        sort!(ogap, [:gap]) #least violated
        ocutid = ogap.id[1:pick]

        global mas = Master(w,sm,1e-4); global dsp = DualSubP(w,1-w,sm,1e-4); MOI.set(mas.m, MOI.LazyConstraintCallback(), Acc_lazy_callback_trunc);
        mvar = [[mas.y[j,t] for j=1:sm.J+sm.K for t=1:2]; [mas.uij[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]]; [mas.ujk[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]]; [mas.uij[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]]; [mas.ujk[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]];[mas.ukl[k,l,m] for k=1:sm.K for l=1:sm.L for m=1:sm.Mkl[k,l]]]

        for o in ocutid
            push!(cutw, [sm.insfile[end-11:end-4], w, dfo.w[o], o])
            @constraint(mas.m, dot((objc.*trunc.(w*dfo.ppi[o][1][1+sum(len[1:6]):end]+(1-w)*dfo.ppi[o][2][1+sum(len[1:6]):end]; digits=dg)), mvar)  <= mas.θ)
        end

        cutime = round(time()-t0; digits=2)
        tCPU = @CPUelapsed optimize!(mas.m)
        if termination_status(mas.m) in [MOI.OPTIMAL,MOI.LOCALLY_SOLVED]
            yt_1,yt_2 = getobjdual(mas.m,dsp.m,sm) 
            val = objective_value(mas.m) #w*yt_1 + (1-w)*yt_2
            
            if (val + buf < lb) && pt[1][1]< yt_1 < pt[2][1] && pt[2][2] < yt_2 < pt[1][2]#round.([yt_1,yt_2]) ∉ YN && val + 10 < lb #&& f1min < yt_1 && yt_2 < f2max#[yt_1,yt_2] ∉ YN #
                println(yt_1," ",yt_2)
                insert = "Y"
                push!(YN, [yt_1,yt_2]); push!(L,[pt[1],[yt_1,yt_2]],[[yt_1,yt_2],pt[2]])
            end 
            # ocutid = oid
            push!(tb, (sm.insfile[end-11:end-4], w, yt_1, yt_2, length(ocutid),sum(tmpnf), sum(tmpno), cutime, sum(subtime), tCPU, insert)) # push!(tb, (sm.insfile[end-11:end-4],, w, lb,val, objective_value(mas.m),yt_1,yt_2, sum(tmpnf), sum(tmpno), sum(tmpit),sum(subtime), tCPU))
            push!(tb2, [pt[1][1],pt[2][1],w,lb, w*yt_1 + (1-w)*yt_2 ,yt_1,yt_2,objective_value(mas.m),insert]) 
        end
        # println("# ocuts: ", size(dfo,1))   
        @label Nextiter
    end
    return YN,tb,tb2,cutw 
end
################## Comparison with the same (fixed) weight #########################

function CBD_fixedW(sm,wdg,rort,buf,f1min,f2max,f1max,f2min,fixedW)
    # global fcuts = []; global ocuts = [] ;
    tb = DataFrame(fname = [], w = [], obj1 = [], obj2 = [], nfcut = [], nocut = [], subCPU = [] ,MPCPU=[]);
    tb2 = DataFrame(w = [], lb = [], val = [], obj1 = [] ,obj2=[], mobj=[]);

    global subtime,tmpnf,tmpno,tmpit = [.0],[0],[0],[0]; global YN = [[f1min,f2max],[f1max,f2min]];
    λ1 = abs(f2max-f2min);  λ2 = abs(f1max-f1min); ratio = λ2/λ1; 
    if rort == "r"
        global w = round(1/(1+ratio); digits=wdg)
    else
        global w = trunc(1/(1+ratio); digits=wdg)
    end
    # println("w = ", w)
    global mas = Master(w,sm,1e-4); global dsp = DualSubP(w,1-w,sm,1e-4); MOI.set(mas.m, MOI.LazyConstraintCallback(), Classic_lazy_callback);
    tCPU = @CPUelapsed optimize!(mas.m); yt_1,yt_2 = getobjdual(mas.m,dsp.m,sm); 
    push!(YN, [yt_1,yt_2]); 
    push!(tb, (sm.insfile[end-11:end-4], w, yt_1, yt_2, sum(tmpnf), sum(tmpno), sum(subtime), tCPU))
    push!(tb2, [w,0,0,yt_1,yt_2,objective_value(mas.m)])

    for w in fixedW
        global (subtime,tmpnf,tmpno,tmpit) = [.0],[0],[0],[0]; 
        # println("w = ", w)
        if w ∈ tb.w || round(w;digits=3)>=.99 || trunc(w;digits=2) <= .01 || w == NaN  
            @goto Nextiter 
        end

        global mas = Master(w,sm,1e-4); global dsp = DualSubP(w,1-w,sm,1e-4); 
        MOI.set(mas.m, MOI.LazyConstraintCallback(), Classic_lazy_callback);
        tCPU = @CPUelapsed optimize!(mas.m)
        if termination_status(mas.m) in [MOI.OPTIMAL,MOI.LOCALLY_SOLVED]
            yt_1,yt_2 = getobjdual(mas.m,dsp.m,sm) 
            println(yt_1," ",yt_2)
            push!(YN, [yt_1,yt_2]);
            push!(tb, (sm.insfile[end-11:end-4], w, yt_1, yt_2, sum(tmpnf), sum(tmpno), sum(subtime), tCPU)) # push!(tb, (sm.insfile[end-11:end-4],, w, lb,val, objective_value(mas.m),yt_1,yt_2, sum(tmpnf), sum(tmpno), sum(tmpit),sum(subtime), tCPU))
            push!(tb2, [w,0, w*yt_1 + (1-w)*yt_2 ,yt_1,yt_2,objective_value(mas.m)]) 
        end    
        @label Nextiter
    end

    return YN,tb,tb2 
end
function CBD_trunc_fixedW(sm,wdg,rort,buf,f1min,f2max,f1max,f2min,fixedW)
    # global fcuts = []; global ocuts = [] ;
    tb = DataFrame(fname = [], w = [], obj1 = [], obj2 = [], nfcut = [], nocut = [], subCPU = [] ,MPCPU=[]);
    tb2 = DataFrame(w = [], lb = [], val = [], obj1 = [] ,obj2=[], mobj=[]);

    global subtime,tmpnf,tmpno,tmpit = [.0],[0],[0],[0]; global YN = [[f1min,f2max],[f1max,f2min]];
    λ1 = abs(f2max-f2min);  λ2 = abs(f1max-f1min); ratio = λ2/λ1; 
    if rort == "r"
        global w = round(1/(1+ratio); digits=wdg)
    else
        global w = trunc(1/(1+ratio); digits=wdg)
    end    
    # println("w = ", w)
    global mas = Master(w,sm,1e-4); global dsp = DualSubP(w,1-w,sm,1e-4); MOI.set(mas.m, MOI.LazyConstraintCallback(), Classic_lazy_callback_trunc);
    tCPU = @CPUelapsed optimize!(mas.m); yt_1,yt_2 = getobjdual(mas.m,dsp.m,sm); 
    lb = w*f1min + (1-w)*f2max;  val = w*yt_1 + (1-w)*yt_2
    push!(YN, [yt_1,yt_2]); 
    push!(tb2, [w,lb,val,yt_1,yt_2,objective_value(mas.m)])
    push!(tb, (sm.insfile[end-11:end-4], w, yt_1, yt_2, sum(tmpnf), sum(tmpno), sum(subtime), tCPU))
    for w in fixedW
        global (subtime,tmpnf,tmpno,tmpit) = [.0],[0],[0],[0]; 
        # println("w = ", w)
        if λ1 ==0 || w ∈ tb.w || round(w;digits=3)>=.991 || trunc(w;digits=2) <= .01 || w == NaN  
            @goto Nextiter 
        end

        global mas = Master(w,sm,1e-4); global dsp = DualSubP(w,1-w,sm,1e-4); 
        MOI.set(mas.m, MOI.LazyConstraintCallback(), Classic_lazy_callback_trunc);
        tCPU = @CPUelapsed optimize!(mas.m)
        if termination_status(mas.m) in [MOI.OPTIMAL,MOI.LOCALLY_SOLVED]
            yt_1,yt_2 = getobjdual(mas.m,dsp.m,sm) 
            println(yt_1," ",yt_2)
            push!(YN, [yt_1,yt_2]);
            push!(tb2, [w,lb, w*yt_1 + (1-w)*yt_2 ,yt_1,yt_2,objective_value(mas.m)]) 
            push!(tb, (sm.insfile[end-11:end-4], w, yt_1, yt_2, sum(tmpnf), sum(tmpno), sum(subtime), tCPU)) # push!(tb, (sm.insfile[end-11:end-4],, w, lb,val, objective_value(mas.m),yt_1,yt_2, sum(tmpnf), sum(tmpno), sum(tmpit),sum(subtime), tCPU))
        end    
        @label Nextiter
    end
    return YN,tb,tb2 
end
function AccBD_fixedW(sm,pcut,wdg,rort,dg,buf,f1min,f2max,f1max,f2min,fixedW)
    tb = DataFrame(fname = [], w = [], obj1 = [], obj2 = [], preOcut = [], nfcut = [], nocut = [], cutCPU =[], subCPU = [] ,MPCPU=[]); 
    tb2 = DataFrame(w = [], lb = [], val = [], obj1 = [] ,obj2=[], mobj=[]);
    global dfo = DataFrame(w=[], ppi=[], α5=[],α6=[]); global dff = DataFrame(w=[],cuts=[])
    global (subtime,tmpnf,tmpno) = ([.0],[0],[0]); global YN = [[f1min,f2max],[f1max,f2min]]; 
    λ1 = abs(f2max-f2min);  λ2 = abs(f1max-f1min); ratio = λ2/λ1; 
    if rort == "r"
        global w = round(1/(1+ratio); digits=wdg)
    else
        global w = trunc(1/(1+ratio); digits=wdg)
    end
    # println("w = ", w)
    global mas = Master(w,sm,1e-4); global dsp = DualSubP(w,1-w,sm,1e-4); MOI.set(mas.m, MOI.LazyConstraintCallback(), Acc_lazy_callback);
    tCPU = @CPUelapsed optimize!(mas.m); yt_1,yt_2 = getobjdual(mas.m,dsp.m,sm); 
    
    push!(YN, [yt_1,yt_2])# println(yt_1," ",yt_2)
    push!(tb, (sm.insfile[end-11:end-4], w, yt_1, yt_2, 0, sum(tmpnf), sum(tmpno), 0, sum(subtime), tCPU))
    push!(tb2, [w,0,0,yt_1,yt_2,objective_value(mas.m)])
    while w in fixedW
        subtime,tmpnf,tmpno = [.0],[0],[0]; 
        if w ∈ tb.w || round(w;digits=3)>=.991 || trunc(w;digits=2) <= .01 || w == NaN  
            @goto Nextiter 
        end
        # Add only violated cuts  
        oid = []; ogap = DataFrame(id=[], gap=[]); #fid = []
        mpobj = value(mas.θ); 
        vy = value.(mas.m[:y]); vuij = value.(mas.m[:uij]); vujk = value.(mas.m[:ujk]); vukl = value.(mas.m[:ukl]);
        masv = [[value.(mas.m[:y])[j,t] for t=1:2 for j=1:sm.J];[value.(mas.m[:y])[k,t] for t=1:2 for k=sm.J+1:sm.J+sm.K];
            [value.(mas.m[:uij])[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]];
            [value.(mas.m[:ujk])[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]];
            [value.(mas.m[:uij])[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]];
            [value.(mas.m[:ujk])[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]];
            [value.(mas.m[:ukl])[k,l,m] for k=1:sm.K for l=1:sm.L for m=1:sm.Mkl[k,l]]]

        t0 = time()       
        global mas = Master(w,sm,1e-4); global dsp = DualSubP(w,1-w,sm,1e-4); MOI.set(mas.m, MOI.LazyConstraintCallback(), Acc_lazy_callback);
        mvar = [[mas.y[j,t] for j=1:sm.J+sm.K for t=1:2]; [mas.uij[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]]; [mas.ujk[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]]; [mas.uij[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]]; [mas.ujk[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]];[mas.ukl[k,l,m] for k=1:sm.K for l=1:sm.L for m=1:sm.Mkl[k,l]]]

        for o=1:size(dfo,1)
            ov = dot(sm.demand,dfo.α5[o])-dot(sm.cas,dfo.α6[o])+dot((objc.*trunc.(w*dfo.ppi[o][1][1+sum(len[1:6]):end]+(1-w)*dfo.ppi[o][2][1+sum(len[1:6]):end]; digits=dg)), masv)
            if ov > mpobj && all(i->i>=0,(w*dfo.ppi[o][1]+(1-w)*dfo.ppi[o][2])) == true            
                push!(oid,o); push!(ogap, (o, abs(ov-mpobj)) )
            end
        end

        # select from sorted violated cuts
        pick = min(pcut,length(oid))
        # sort!(ogap, [:gap], rev=true)
        sort!(ogap, [:gap])
        ocutid = ogap.id[1:pick]
        for o in ocutid
            @constraint(mas.m, dot((objc.*trunc.(w*dfo.ppi[o][1][1+sum(len[1:6]):end]+(1-w)*dfo.ppi[o][2][1+sum(len[1:6]):end]; digits=dg)), mvar)  <= mas.θ)
        end
        
        cutime = round(time()-t0; digits=2)
        tCPU = @CPUelapsed optimize!(mas.m)
        if termination_status(mas.m) in [MOI.OPTIMAL,MOI.LOCALLY_SOLVED]
            yt_1,yt_2 = getobjdual(mas.m,dsp.m,sm) 
            println(yt_1," ",yt_2)
            push!(YN, [yt_1,yt_2]);
            #ocutid = oid
            push!(tb2, [w,0, w*yt_1 + (1-w)*yt_2 ,yt_1,yt_2,objective_value(mas.m)]) 
            push!(tb, (sm.insfile[end-11:end-4], w, yt_1, yt_2, length(ocutid),sum(tmpnf), sum(tmpno), cutime, sum(subtime), tCPU, insert)) 
        end   
        @label Nextiter
    end
    return YN,tb,tb2
end
function AccBD_trunc_fixedW(sm,pcut,wdg,rort,dg,buf,f1min,f2max,f1max,f2min,fixedW)
    tb = DataFrame(fname = [], w = [], obj1 = [], obj2 = [], preOcut = [], nfcut = [], nocut = [], cutCPU =[], subCPU = [] ,MPCPU=[]);   
    tb2 = DataFrame(w = [], lb = [], val = [], obj1 = [] ,obj2=[], mobj=[]);

    global dfo = DataFrame(w=[], ppi=[], α5=[],α6=[]); global dff = DataFrame(w=[],cuts=[])# global (pi1,pi2,fcuts,ocuts) = [],[],[],[] 
    global (subtime,tmpnf,tmpno,tmpit) = ([.0],[0],[0],[0]); global YN = [[f1min,f2max],[f1max,f2min]]; # global ocuts = []; 
    λ1 = abs(f2max-f2min);  λ2 = abs(f1max-f1min); ratio = λ2/λ1; 
    if rort == "r"
        global w = round(1/(1+ratio); digits=wdg)
    else
        global w = trunc(1/(1+ratio); digits=wdg)
    end
    # println("w = ", w);  
    global mas = Master(w,sm,1e-4); global dsp = DualSubP(w,1-w,sm,1e-4); MOI.set(mas.m, MOI.LazyConstraintCallback(), Acc_lazy_callback_trunc);
    tCPU = @CPUelapsed optimize!(mas.m); yt_1,yt_2 = getobjdual(mas.m,dsp.m,sm); 
    push!(YN, [yt_1,yt_2])# println(yt_1," ",yt_2)
    push!(tb2, [w,0,0,yt_1,yt_2,objective_value(mas.m)])
    push!(tb, (sm.insfile[end-11:end-4], w, yt_1, yt_2, 0, sum(tmpnf), sum(tmpno), 0, sum(subtime), tCPU))
    for w in fixedW
        subtime,tmpnf,tmpno,tmpit = [.0],[0],[0],[0]; 
        # println(" w = ", w) 
        if w ∈ tb.w || round(w;digits=3)>=.99 || trunc(w;digits=2) <= .01 || w == NaN
            @goto Nextiter
        end

        # Add only violated cuts
        oid = []; ogap = DataFrame(id=[], gap=[])    #fid = []
        mpobj = value(mas.θ);
        masv = [[value.(mas.m[:y])[j,t] for t=1:2 for j=1:sm.J];[value.(mas.m[:y])[k,t] for t=1:2 for k=sm.J+1:sm.J+sm.K];
            [value.(mas.m[:uij])[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]];
            [value.(mas.m[:ujk])[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]];
            [value.(mas.m[:uij])[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]];
            [value.(mas.m[:ujk])[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]];
            [value.(mas.m[:ukl])[k,l,m] for k=1:sm.K for l=1:sm.L for m=1:sm.Mkl[k,l]]]

        t0 = time()
        
        global mas = Master(w,sm,1e-4); global dsp = DualSubP(w,1-w,sm,1e-4); MOI.set(mas.m, MOI.LazyConstraintCallback(), Acc_lazy_callback_trunc);
        mvar = [[mas.y[j,t] for j=1:sm.J+sm.K for t=1:2]; [mas.uij[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]]; [mas.ujk[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]]; [mas.uij[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]]; [mas.ujk[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]];[mas.ukl[k,l,m] for k=1:sm.K for l=1:sm.L for m=1:sm.Mkl[k,l]]]

        for o=1:size(dfo,1)
            ov = dot(sm.demand,dfo.α5[o]) - dot(sm.cas,dfo.α6[o]) + dot((objc.*trunc.(w*dfo.ppi[o][1][1+sum(len[1:6]):end]+(1-w)*dfo.ppi[o][2][1+sum(len[1:6]):end]; digits=dg)), masv)
            if ov > mpobj && all(i->i>=0,(w*dfo.ppi[o][1]+(1-w)*dfo.ppi[o][2])) == true
                push!(oid,o); push!(ogap, (o, abs(ov-mpobj)) )
            end
        end
        #select from sorted violated cuts
        pick = min(pcut,length(oid))
        # sort!(ogap, [:gap], rev=true)
        sort!(ogap, [:gap])
        ocutid = ogap.id[1:pick]
        for o in ocutid
            @constraint(mas.m, dot((objc.*trunc.(w*dfo.ppi[o][1][1+sum(len[1:6]):end]+(1-w)*dfo.ppi[o][2][1+sum(len[1:6]):end]; digits=dg)), mvar)  <= mas.θ)
        end

        cutime = round(time()-t0; digits=2)
        tCPU = @CPUelapsed optimize!(mas.m)
        if termination_status(mas.m) in [MOI.OPTIMAL,MOI.LOCALLY_SOLVED]

            yt_1,yt_2 = getobjdual(mas.m,dsp.m,sm) 
            println(yt_1," ",yt_2)
            push!(YN, [yt_1,yt_2]); 
            # ocutid = oid
            push!(tb2, [w,0, w*yt_1 + (1-w)*yt_2 ,yt_1,yt_2,objective_value(mas.m)]) 
            push!(tb, (sm.insfile[end-11:end-4], w, yt_1, yt_2, length(ocutid),sum(tmpnf), sum(tmpno), cutime, sum(subtime), tCPU)) # push!(tb, (sm.insfile[end-11:end-4],, w, lb,val, objective_value(mas.m),yt_1,yt_2, sum(tmpnf), sum(tmpno), sum(tmpit),sum(subtime), tCPU))
        end
        @label Nextiter
    end
    return YN,tb,tb2
end

