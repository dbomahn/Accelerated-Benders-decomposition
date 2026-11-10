# sm  = Data(file);
# file = "../instances/small/test01S1.dat"; tt = GenerateCuts(saveto,file,tol)

include(joinpath(@__DIR__, "paths.jl"))
include(joinpath(@__DIR__, "..", "src", "AcceleratedBD.jl"))
# global fname = file[end-11:end-4]
# ins = readdir(pathh)[1:end]
pathh = INSTANCE_DIR
ins = readdir(pathh)[2:end]
i = 1
sm = Data(joinpath(pathh, ins[i])); pcut= 40;
if sm.insfile[end-11:end-4] == "test01S1"
    rort = "t"; global num = 5; buf = 10^4; wdg = 4; 
elseif sm.insfile[end-11:end-4] == "test01S2"
    rort = "r"; global num = 3; buf = 10^4; wdg = 4;
elseif sm.insfile[end-11:end-4] == "test01S3"
    rort = "r"; global num = 4; buf = 14000; wdg = 4;
elseif sm.insfile[end-11:end-4] == "test01S4"
    rort = "t"; global num = 4; buf = 9000;  wdg = 4;
elseif sm.insfile[end-11:end-4] == "test02S1"
    rort = "r"; global num = 2; buf = 10^4; wdg = 3;
elseif sm.insfile[end-11:end-4] == "test02S2" 
    rort = "r"; global num = 3; buf = 10^4; wdg = 4;
elseif sm.insfile[end-11:end-4] == "test02S3"
    rort = "r"; global num = 4; buf = 10^4; wdg = 4;
elseif sm.insfile[end-11:end-4] == "test02S4" 
    rort = "r"; global num = 4; buf = 10^4; wdg = 4;
elseif sm.insfile[end-11:end-4] == "test03S2" 
    rort = "r"; global num = 4; buf = 10^4; wdg = 4;
else #sm.insfile[end-11:end-4] == "test03S3"
    rort = "r"; global num = 4; buf = 8000; wdg = 4;
end       

dtime = @CPUelapsed vmodel = smSCND_MIP(sm,1e-4); vtime = @CPUelapsed vSolve(vmodel,Inf,4,method=:dicho,verbose=false); res = unique(getY_N(vmodel)); res; v1 = [res[i][1] for i=1:length(res)]; v2 = [res[i][2] for i=1:length(res)]; v = scatter(x=v1, y=v2, name="vopt", mode="markers", market=attr(color="lime"))
plot([v],layout)
# clipboard(dtime)
# cdicho = @CPUelapsed csol,ct,ctt = ClassicBD(sm, wdg, rort, buf); adicho = @CPUelapsed asol,at,att = AcceleratedBD(sm, pcut, wdg, rort, 1, buf); cy = NDfilter(csol); csol1 = [csol[i][1] for i=1:length(csol)]; csol2 = [csol[i][2] for i=1:length(csol)]; t3 = scatter(x=csol1, y=csol2, name="Classic", mode="markers", market=attr(color="red")); ay = NDfilter(asol); asol1= [asol[i][1] for i=1:length(asol)]; asol2 = [asol[i][2] for i=1:length(asol)]; t1 = scatter(x=asol1, y=asol2, name="Acc", mode="markers", marker=attr(color="lime")); println("round ", ct,at);plot([v,t3,t1],layout)
# plot([v,t3,t1],layout)

#########
cdicho = @CPUelapsed csol,ct,ctt = ClassicBD_trunc(sm, wdg, rort, buf); 
cy = NDfilter(csol); csol1 = [cy[i][1] for i=1:length(cy)]; csol2 = [cy[i][2] for i=1:length(cy)]; t3 = scatter(x=csol1, y=csol2, name="Classic", mode="markers", market=attr(color="red")); 
plot([v,t3],layout)
CSV.write(joinpath(RESULTS_DIR, "T1S1_CB.xlsx"), at, append=true, header=true, delim=' '); 

adicho = @CPUelapsed asol,at,att = AcceleratedBD_trunc(sm, pcut, wdg, rort, 1, buf); 
ay = NDfilter(asol); asol1= [ay[i][1] for i=1:length(asol)]; asol2 = [ay[i][2] for i=1:length(ay)]; t1 = scatter(x=asol1, y=asol2, name="Acc", mode="markers", marker=attr(color="lime")); println("truncate ",ct,at); plot([v,t3,t1],layout)
CSV.write(joinpath(RESULTS_DIR, "ABDresults_leastViolated3.xlsx"), at, append=true, header=true, delim=' ')
CSV.write(joinpath(RESULTS_DIR, "ABDtime.xlsx"), DataFrame(ins = sm.insfile[end-11:end-4], cput = round(adicho; digits=1), cut = "Lviolated"), append=true )

plot([v,t3,t1],layout)
i=2
at.w[i]*at.obj1[i] + (1-at.w[i])*at.obj2[i]
ct.w[i]*ct.obj1[i] + (1-ct.w[i])*ct.obj2[i]

for k=1:5
    @show k
    if i in [1,8,9]
        cdicho = @CPUelapsed csol,ct,ctt = ClassicBD_trunc(sm, 4, rort, buf); 
        println(ct) 
        adicho = @CPUelapsed asol,at,att = AcceleratedBD_trunc(sm, pcut, 4, rort, 1, buf); println(at)
        CSV.write(joinpath(RESULTS_DIR, "ABDtime$pcut.xlsx"), DataFrame(ins = sm.insfile[end-11:end-4], cput = round(adicho; digits=1), cut = "Lviolated_$pcut"), append=true )
        CSV.write(joinpath(RESULTS_DIR, "AB$i")*"_"*"$pcut.xlsx", at, append=true, header=true, delim=' ')
        CSV.write(joinpath(RESULTS_DIR, "CBDtime$pcut.xlsx"), DataFrame(ins = sm.insfile[end-11:end-4], cput = round(cdicho; digits=1)), append=true )
        CSV.write(joinpath(RESULTS_DIR, "CB$i")*"_"*"$pcut.xlsx", ct, append=true, header=true, delim=' '); 
    else
        cdicho = @CPUelapsed csol,ct,ctt = ClassicBD(sm, wdg, rort, buf); 
        println(ct)
        adicho = @CPUelapsed asol,at,att = AcceleratedBD(sm, pcut, wdg, rort, 1, buf); 
        println(at)
        CSV.write(joinpath(RESULTS_DIR, "ABDtime$pcut.xlsx"), DataFrame(ins = sm.insfile[end-11:end-4], cput = round(adicho; digits=1), cut = "Lviolated_$pcut"), append=true )
        CSV.write(joinpath(RESULTS_DIR, "AB$i")*"_"*"$pcut.xlsx", at, append=true, header=true, delim=' ')
        CSV.write(joinpath(RESULTS_DIR, "CBDtime$pcut.xlsx"), DataFrame(ins = sm.insfile[end-11:end-4], cput = round(cdicho; digits=1)), append=true )
        CSV.write(joinpath(RESULTS_DIR, "CB$i")*"_"*"$pcut.xlsx", ct, append=true, header=true, delim=' '); 
        # CSV.write(joinpath(RESULTS_DIR, "CutsFrom.xlsx"), ogw, append=true , header=true, delim=' ')
    end
end

for i=1:5
    cdicho = @CPUelapsed csol,ct,ctt = ClassicBD(sm, 4, rort, buf); 
    adicho = @CPUelapsed asol,at,att = AcceleratedBD(sm, pcut, 4, rort, 1, buf); 
    # cdicho = @CPUelapsed csol,ct,ctt = ClassicBD_trunc(sm, 4, rort, buf); 
    # adicho = @CPUelapsed asol,at,att = AcceleratedBD_trunc(sm, pcut, 4, rort, 1, buf); 
    CSV.write(joinpath(RESULTS_DIR, "T2S3.xlsx"), ct, append=true, header=true, delim=' ')
    CSV.write(joinpath(RESULTS_DIR, "T2S3.xlsx"), at, append=true, header=true, delim=' ')
    CSV.write(joinpath(RESULTS_DIR, "BD_sameweight_time.xlsx"), DataFrame(ctime=cdicho, atime=adicho), append=true, header=false, delim=' ')
end
CSV.write(joinpath(RESULTS_DIR, "BD_sameweight_time.xlsx"), DataFrame(ctime=cdicho, atime=adicho), append=true, header=false, delim=' ')
CSV.write(joinpath(RESULTS_DIR, "T1S1.xlsx"), at, append=true, header=true, delim=' ')


# Truncate
cdicho = @CPUelapsed csol,ct,ctt = ClassicBD_trunc(sm, 4, rort, buf); println(ct)
CSV.write(joinpath(RESULTS_DIR, "CBDresults3.xlsx"), ct, append=true, header=true, delim=' '); 
CSV.write(joinpath(RESULTS_DIR, "CBDtime.xlsx"), DataFrame(ins = sm.insfile[end-11:end-4], cput = round(cdicho; digits=1)), append=true )

for pcut in [10,20,30,40,50]#[60,80,100,120]#
    adicho = @CPUelapsed asol,at,att = AcceleratedBD_trunc(sm, pcut, 4, rort, 1, buf); println(at)  
    CSV.write(joinpath(RESULTS_DIR, "ABDresults_leastViolated2.xlsx"), at, append=true, header=true, delim=' ')
    CSV.write(joinpath(RESULTS_DIR, "ABDtime.xlsx"), DataFrame(ins = sm.insfile[end-11:end-4], cput = round(adicho; digits=1), cut = "Lviolated_$pcut"), append=true )
end
# Round
cdicho = @CPUelapsed csol,ct,ctt = ClassicBD(sm, 4, rort, buf); println(ct)
CSV.write(joinpath(RESULTS_DIR, "CBDresults3.xlsx"), ct, append=true, header=true, delim=' ')
CSV.write(joinpath(RESULTS_DIR, "CBDtime.xlsx"), DataFrame(ins = sm.insfile[end-11:end-4], cput = round(cdicho; digits=1)), append=true )

for pcut in [10,20,30,40,50]#[60,70,80,100,120]
    adicho = @CPUelapsed asol,at,att = AcceleratedBD(sm, pcut, 4, rort, 1, buf); #println(at)
    CSV.write(joinpath(RESULTS_DIR, "ABDresults_leastViolated2.xlsx"), at, append=true, header=true, delim=' ')
    CSV.write(joinpath(RESULTS_DIR, "ABDtime.xlsx"), DataFrame(ins = sm.insfile[end-11:end-4], cput = round(adicho; digits=1), cut = "Lviolated_$pcut"), append=true )
end1


# for f=1:size(dff,1) #length(fcuts) #
#     fv = dot(sm.demand,dff.cuts[f].α5)-dot(sm.cas,dff.cuts[f].α6) -
#     sum(sm.cap[j]*dff.cuts[f].α7[j,t]*vy[j,t] for j=1:sm.J for t=1:2) -
#     sum(sm.cad[k]*dff.cuts[f].α8[k,t]*vy[k+sm.J,t] for k=1:sm.K for t=1:2) +
#     sum(sm.Lcapasp[i][j][m]*dff.cuts[f].α9[i,j,m]*vuij[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]) +
#     sum(sm.Lcapapd[j][k][m]*dff.cuts[f].α10[j,k,m]*vujk[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k])-
#     sum(sm.bigM*vuij[i,j,m]*dff.cuts[f].α12[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]) -
#     sum(sm.bigM*vujk[j,k,m]*dff.cuts[f].α13[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]) -
#     sum(sm.bigM*vukl[k,l,m]*dff.cuts[f].α14[k,l,m] for k=1:sm.K for l=1:sm.L for m=1:sm.Mkl[k,l])
#     # fv = dot(sm.demand,dff.cuts[f].α5)-dot(sm.cas,dff.cuts[f].α6)- 
#     # sum(sm.cap[j]*trunc.(dff.cuts[f].α7[j,t]; digits=dg)*vy[j,t] for j=1:sm.J for t=1:2)- 
#     # sum(sm.cad[k]*trunc.(dff.cuts[f].α8[k,t]; digits=dg)*vy[k+sm.J,t] for k=1:sm.K for t=1:2) + 
#     # sum(sm.Lcapasp[i][j][m]*trunc.(dff.cuts[f].α9[i,j,m]; digits=dg)*vuij[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]) +
#     # sum(sm.Lcapapd[j][k][m]*trunc.(dff.cuts[f].α10[j,k,m]; digits=dg)*vujk[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k])-
#     # sum(sm.bigM*vuij[i,j,m]*trunc.(dff.cuts[f].α12[i,j,m]; digits=dg) for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]) -
#     # sum(sm.bigM*vujk[j,k,m]*trunc.(dff.cuts[f].α13[j,k,m]; digits=dg) for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]) -
#     # sum(sm.bigM*vukl[k,l,m]*trunc.(dff.cuts[f].α14[k,l,m]) for k=1:sm.K for l=1:sm.L for m=1:sm.Mkl[k,l])

#     if fv > 0
#         # println("feasibility cut = ", fv)                
#         push!(fid,f) # println(fid)
#     end
# end
# ov = dot(sm.demand,dfo.α5[o])-dot(sm.cas,dfo.α6[o])+dot((objc.*round.(w*dfo.ppi[o][1][1+sum(len[1:6]):end]+(1-w)*dfo.ppi[o][2][1+sum(len[1:6]):end]; digits=dg)), masv)
# ov = dot(sm.demand,dfo.α5[o])-dot(sm.cas,dfo.α6[o]) -
#     dot([sm.cap[j] for t=1:2 for j=1:sm.J].*(w*dfo.ppi[o][1][1+sum(len[1:6]):sum(len[1:7])]+(1-w)*dfo.ppi[o][2][1+sum(len[1:6]):sum(len[1:7])]),[vy[j,t] for t=1:2 for j=1:sm.J]) -
#     dot([sm.cad[k] for t=1:2 for k=1:sm.K].*(w*dfo.ppi[o][1][1+sum(len[1:7]):sum(len[1:8])]+(1-w)*dfo.ppi[o][2][1+sum(len[1:7]):sum(len[1:8])]),[vy[k+sm.J,t] for t=1:2 for k=1:sm.K])+ 
#     dot([sm.Lcapasp[i][j][m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]].*(w*dfo.ppi[o][1][1+sum(len[1:8]):sum(len[1:9])]+(1-w)*dfo.ppi[o][2][1+sum(len[1:8]):sum(len[1:9])]),[vuij[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]])+ 
#     dot([sm.Lcapapd[j][k][m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]].*(w*dfo.ppi[o][1][1+sum(len[1:9]):sum(len[1:10])]+(1-w)*dfo.ppi[o][2][1+sum(len[1:9]):sum(len[1:10])]),[vujk[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]])-
#     dot([sm.bigM for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]].*(w*dfo.ppi[o][1][1+sum(len[1:10]):sum(len[1:11])]+(1-w)*dfo.ppi[o][2][1+sum(len[1:10]):sum(len[1:11])]),[vuij[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]])-
#     dot([sm.bigM for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]].*(w*dfo.ppi[o][1][1+sum(len[1:11]):sum(len[1:12])]+(1-w)*dfo.ppi[o][2][1+sum(len[1:11]):sum(len[1:12])]), [vujk[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]])- 
#     dot([sm.bigM for k=1:sm.K for l=1:sm.L for m=1:sm.Mkl[k,l]].*(w*dfo.ppi[o][1][1+sum(len[1:12]):end]+(1-w)*dfo.ppi[o][2][1+sum(len[1:12]):end]), [vukl[k,l,m] for k=1:sm.K for l=1:sm.L for m=1:sm.Mkl[k,l]])
atb = []; 
for num in [10,20,30,40]#[20,40,60,80]# 
    # adicho = @CPUelapsed asol,at = AcceleratedBD(sm,4,"r",num,num);push!(atb,adicho)
    adicho = @CPUelapsed asol,at = AcceleratedBD(sm,4,"r",2*num,0);push!(atb,adicho)
    # adicho = @CPUelapsed asol,at = AcceleratedBD(sm,4,"r",0,2*num);push!(atb,adicho)

    # adicho = @CPUelapsed asol,at = AcceleratedBD(sm,4,"t",num,num);push!(atb,adicho)
    # adicho = @CPUelapsed asol,at = AcceleratedBD(sm,4,"t",2*num,0);push!(atb,adicho)
    # adicho = @CPUelapsed asol,at = AcceleratedBD(sm,4,"t",0,2*num);push!(atb,adicho)

    # CSV.write(joinpath(RESULTS_DIR, "Benders_result.csv"), at, append=true, header=true, delim=' ')
end
atb

###
for k=1:10
    cdicho = @CPUelapsed cy,ctb= ClassicBD(ins[k])
    CSV.write(joinpath(RESULTS_DIR, "CBdicho_Iter.csv"), ctb, append=true, header=true, delim=' ')

    cb2 = DataFrame(fname ="$fname", fcut = sum(ctb.fcut), ocut = sum(ctb.ocut), iter = sum(ctb.iter), subT = round(sum(ctb.subCPU);digits=2), BendT = round(sum(ctb.MPCPU);digits=2),DichoT=round(cdicho;digits=2))
    CSV.write(joinpath(RESULTS_DIR, "ClassicBD_tb.csv"), cb2, append=true, header=false, delim=' ')

    ctable = zeros(length(cy),2)
    for i=1:length(cy)
        for j=1:2
            ctable[i,:] = cy[i]
        end
    end
    CSV.write(CLASSIC_NDP_DIR * "/"*fname*"Y.log",DataFrame(ctable, :auto),append=false, header=false, delim=' ' )
end
for j=5:6
    adicho = @CPUelapsed ay,atb = AcceleratedBD(ins[j])
    CSV.write(joinpath(RESULTS_DIR, "ABdicho_Iter.csv"), atb, append=true, header=true, delim=' ')

    ab2 = DataFrame(fname ="$fname", fcut = sum(atb.fcut), ocut = sum(atb.ocut), iter = sum(atb.iter), subT = round(sum(atb.subCPU);digits=2), BendT = round(sum(atb.MPCPU);digits=2),DichoT=round(adicho;digits=2))
    CSV.write(joinpath(RESULTS_DIR, "AccBD_tb.csv"), ab2, append=true, header=false, delim=' ')

    otable = zeros(length(ay),2)
    for i=1:length(ay)
        for j=1:2
            otable[i,:] = ay[i]
        end
    end
    CSV.write(ACCELERATED_NDP_DIR * "/"*fname*"Y.log",DataFrame(otable, :auto),append=false, header=false, delim=' ' )
end

# global mas = Master(w,sm,1e-4); global dsp = DualSubP(w,1-w,sm,1e-4); MOI.set(mas.m, MOI.LazyConstraintCallback(), Acc_lazy_callback);
# mvar = [[mas.y[j,t] for j=1:sm.J+sm.K for t=1:2]; [mas.uij[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]]; [mas.ujk[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]]; [mas.uij[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]]; [mas.ujk[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]];[mas.ukl[k,l,m] for k=1:sm.K for l=1:sm.L for m=1:sm.Mkl[k,l]]]

# @constraints(mas.m, begin 
        # [f in fid], 
        # dot(sm.demand,dff.cuts[f].α5)-dot(sm.cas,dff.cuts[f].α6)- 
        # sum(sm.cap[j]*trunc.(dff.cuts[f].α7[j,t]; digits=dg)*mas.y[j,t] for j=1:sm.J for t=1:2)- 
        # sum(sm.cad[k]*trunc.(dff.cuts[f].α8[k,t]; digits=dg)*mas.y[k+sm.J,t] for k=1:sm.K for t=1:2) + 
        # sum(sm.Lcapasp[i][j][m]*trunc.(dff.cuts[f].α9[i,j,m]; digits=dg)*mas.uij[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]) +
        # sum(sm.Lcapapd[j][k][m]*trunc.(dff.cuts[f].α10[j,k,m]; digits=dg)*mas.ujk[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k])-
        # sum(sm.bigM*mas.uij[i,j,m]*trunc.(dff.cuts[f].α12[i,j,m]; digits=dg) for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]) -
        # sum(sm.bigM*mas.ujk[j,k,m]*trunc.(dff.cuts[f].α13[j,k,m]; digits=dg) for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]) -
        # sum(sm.bigM*mas.ukl[k,l,m]*trunc.(dff.cuts[f].α14[k,l,m]; digits=dg) for k=1:sm.K for l=1:sm.L for m=1:sm.Mkl[k,l]) <= 0

        # dot(sm.demand,dff.cuts[f].α5)-dot(sm.cas,dff.cuts[f].α6)-
        # sum(sm.cap[j]*dff.cuts[f].α7[j,t]*mas.y[j,t] for j=1:sm.J for t=1:2) -
        # sum(sm.cad[k]*dff.cuts[f].α8[k,t]*mas.y[k+sm.J,t] for k=1:sm.K for t=1:2) +
        # sum(sm.Lcapasp[i][j][m]*dff.cuts[f].α9[i,j,m]*mas.uij[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]) +
        # sum(sm.Lcapapd[j][k][m]*dff.cuts[f].α10[j,k,m]*mas.ujk[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]) -
        # sum(sm.bigM*dff.cuts[f].α12[i,j,m]*mas.uij[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]) -
        # sum(sm.bigM*dff.cuts[f].α13[j,k,m]*mas.ujk[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]) -
        # sum(sm.bigM*dff.cuts[f].α14[k,l,m]*mas.ukl[k,l,m] for k=1:sm.K for l=1:sm.L for m=1:sm.Mkl[k,l]) <= 0

        # [o in oid],
        # dot(sm.demand,dfo.α5[o])-dot(sm.cas,dfo.α6[o]) -
        # dot([sm.cap[j] for t=1:2 for j=1:sm.J].*trunc.(w*dfo.ppi[o][1][1+sum(len[1:6]):sum(len[1:7])]+(1-w)*dfo.ppi[o][2][1+sum(len[1:6]):sum(len[1:7])]; digits=dg), mas.y[j,t] for j=1:sm.J for t=1:2) -
        # dot([sm.cad[k] for t=1:2 for k=1:sm.K].*trunc.(w*dfo.ppi[o][1][1+sum(len[1:7]):sum(len[1:8])]+(1-w)*dfo.ppi[o][2][1+sum(len[1:7]):sum(len[1:8])]; digits=dg), mas.y[k+sm.J,t] for k=1:sm.K for t=1:2)+ 
        # dot([sm.Lcapasp[i][j][m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]].*trunc.(w*dfo.ppi[o][1][1+sum(len[1:8]):sum(len[1:9])]+(1-w)*dfo.ppi[o][2][1+sum(len[1:8]):sum(len[1:9])]; digits=dg), mas.uij[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j])+ 
        # dot([sm.Lcapapd[j][k][m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]].*trunc.(w*dfo.ppi[o][1][1+sum(len[1:9]):sum(len[1:10])]+(1-w)*dfo.ppi[o][2][1+sum(len[1:9]):sum(len[1:10])]; digits=dg), mas.ujk[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k])-
        # dot([sm.bigM for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]].*trunc.(w*dfo.ppi[o][1][1+sum(len[1:10]):sum(len[1:11])]+(1-w)*dfo.ppi[o][2][1+sum(len[1:10]):sum(len[1:11])]; digits=dg), mas.uij[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j])-
        # dot([sm.bigM for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]].*trunc.(w*dfo.ppi[o][1][1+sum(len[1:11]):sum(len[1:12])]+(1-w)*dfo.ppi[o][2][1+sum(len[1:11]):sum(len[1:12])]; digits=dg), mas.ujk[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k])- 
        # dot([sm.bigM for k=1:sm.K for l=1:sm.L for m=1:sm.Mkl[k,l]].*trunc.(w*dfo.ppi[o][1][1+sum(len[1:12]):end]+(1-w)*dfo.ppi[o][2][1+sum(len[1:12]):end]; digits=dg), mas.ukl[k,l,m] for k=1:sm.K for l=1:sm.L for m=1:sm.Mkl[k,l]) <= mas.θ                
        
        # dot(sm.demand,dfo.α5[o])-dot(sm.cas,dfo.α6[o]) -
        # dot([sm.cap[j] for t=1:2 for j=1:sm.J].*(w*dfo.ppi[o][1][1+sum(len[1:6]):sum(len[1:7])]+(1-w)*dfo.ppi[o][2][1+sum(len[1:6]):sum(len[1:7])]), mas.y[j,t] for j=1:sm.J for t=1:2) -
        # dot([sm.cad[k] for t=1:2 for k=1:sm.K].*(w*dfo.ppi[o][1][1+sum(len[1:7]):sum(len[1:8])]+(1-w)*dfo.ppi[o][2][1+sum(len[1:7]):sum(len[1:8])]), mas.y[k+sm.J,t] for k=1:sm.K for t=1:2)+ 
        # dot([sm.Lcapasp[i][j][m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]].*(w*dfo.ppi[o][1][1+sum(len[1:8]):sum(len[1:9])]+(1-w)*dfo.ppi[o][2][1+sum(len[1:8]):sum(len[1:9])]), mas.uij[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j])+ 
        # dot([sm.Lcapapd[j][k][m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]].*(w*dfo.ppi[o][1][1+sum(len[1:9]):sum(len[1:10])]+(1-w)*dfo.ppi[o][2][1+sum(len[1:9]):sum(len[1:10])]), mas.ujk[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k])-
        # dot([sm.bigM for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]].*(w*dfo.ppi[o][1][1+sum(len[1:10]):sum(len[1:11])]+(1-w)*dfo.ppi[o][2][1+sum(len[1:10]):sum(len[1:11])]), mas.uij[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j])-
        # dot([sm.bigM for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]].*(w*dfo.ppi[o][1][1+sum(len[1:11]):sum(len[1:12])]+(1-w)*dfo.ppi[o][2][1+sum(len[1:11]):sum(len[1:12])]), mas.ujk[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k])- 
        # dot([sm.bigM for k=1:sm.K for l=1:sm.L for m=1:sm.Mkl[k,l]].*(w*dfo.ppi[o][1][1+sum(len[1:12]):end]+(1-w)*dfo.ppi[o][2][1+sum(len[1:12]):end]), mas.ukl[k,l,m] for k=1:sm.K for l=1:sm.L for m=1:sm.Mkl[k,l]) <= mas.θ
        
        # dot(sm.demand,dfo.α5[o])-dot(sm.cas,dfo.α6[o]) +
        # dot((w*dfo.ppi[o][1][1:numcon]+(1-w)*dfo.ppi[o][2][1:numcon]),all_variables(mas.m)[1:end-1]) <= mas.θ
        
# end)


# ab2 = DataFrame(fname ="$fname", fcut = sum(at.fcut), ocut = sum(at.ocut), 
#     iter = sum(at.iter), subT = round(sum(at.subCPU);digits=2), 
#      BendT = round(sum(at.totalCPU);digits=2),DichoT=round(adicho;digits=2))
# CSV.write(joinpath(RESULTS_DIR, "AccBD_tb.csv"), ab2, append=true, header=false, delim=' ')

# otable = zeros(length(ay),2)
# for i=1:length(ay)
#     for j=1:2
#         otable[i,:] = ay[i]
#     end
# end
# CSV.write(ACCELERATED_NDP_DIR * "/"*fname*"Y.log",DataFrame(otable, :auto),append=false, header=false, delim=' ' )

# Add cuts after sampling  
# fcutid = sample(1:length(fcuts), min(length(fcuts),sm.I*5), replace = false)
# ocutid = sample(1:length(pi1), min(length(pi1),sm.I*5), replace = false) #sm.I*5 
# fcutid = sample(1:length(fcuts), min(length(fcuts),fnum), replace = false)
# ocutid = sample(1:length(pi1), min(length(pi1),onum), replace = false) #sm.I*5 
# ocutid = sample(1:length(ocuts), min(length(ocuts),onum), replace = false) #sm.I*5 
# Add cuts from pre-generated pool
# ocutid,fcutid = CutSelection_closweight(w,dfo,dff,onum,fnum) #o1id,o2id,f1id,f2id; ocutid = [o1id;o2id]; fcutid = [f1id;f2id];
# global mas = Master(w,sm,1e-4); global dsp = DualSubP(w,1-w,sm,1e-4);
# MOI.set(mas.m, MOI.LazyConstraintCallback(), Acc_lazy_callback);
# @constraints(mas.m, begin
#     [f in fcutid],# [f = 1:length(fcuts)],#
#     dot(sm.demand,dff.cuts[f].α5)-dot(sm.cas,dff.cuts[f].α6)-
#     sum(sm.cap[j]*dff.cuts[f].α7[j,t]*mas.y[j,t] for j=1:sm.J for t=1:2) -
#     sum(sm.cad[k]*dff.cuts[f].α8[k,t]*mas.y[k+sm.J,t] for k=1:sm.K for t=1:2) +
#     sum(sm.Lcapasp[i][j][m]*dff.cuts[f].α9[i,j,m]*mas.uij[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]) +
#     sum(sm.Lcapapd[j][k][m]*dff.cuts[f].α10[j,k,m]*mas.ujk[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]) -
#     sum(sm.bigM*mas.uij[i,j,m]*dff.cuts[f].α12[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]) -
#     sum(sm.bigM*mas.ujk[j,k,m]*dff.cuts[f].α13[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]) -
#     sum(sm.bigM*mas.ukl[k,l,m]*dff.cuts[f].α14[k,l,m] for k=1:sm.K for l=1:sm.L for m=1:sm.Mkl[k,l]) <= 0
#     #  dot(sm.demand,dff.cuts[f].α5)-dot(sm.cas,dff.cuts[f].α6)-
#     #  sum(sm.cap[j]*dff.cuts[f].α7[j,t]*mas.y[j,t] for j=1:sm.J for t=1:2) -
#     #  sum(sm.cad[k]*dff.cuts[f].α8[k,t]*mas.y[k+sm.J,t] for k=1:sm.K for t=1:2) +
#     #  sum(sm.Lcapasp[i][j][m]*dff.cuts[f].α9[i,j,m]*mas.uij[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]) +
#     #  sum(sm.Lcapapd[j][k][m]*dff.cuts[f].α10[j,k,m]*mas.ujk[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]) -
#     #  sum(sm.bigM*mas.uij[i,j,m]*dff.cuts[f].α12[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]) -
#     #  sum(sm.bigM*mas.ujk[j,k,m]*dff.cuts[f].α13[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]) -
#     #  sum(sm.bigM*mas.ukl[k,l,m]*dff.cuts[f].α14[k,l,m] for k=1:sm.K for l=1:sm.L for m=1:sm.Mkl[k,l]) <= 0
#     [o in ocutid], # [o=1:length(pi1)],#
#         dot(sm.demand,ocuts[o].α5)-dot(sm.cas,ocuts[o].α6)-
#         sum(sm.cap[j]*ocuts[o].α7[j,t]*mas.y[j,t] for j=1:sm.J for t=1:2) -
#         sum(sm.cad[k]*ocuts[o].α8[k,t]*mas.y[k+sm.J,t] for k=1:sm.K for t=1:2) +
#         sum(sm.Lcapasp[i][j][m]*ocuts[o].α9[i,j,m]*mas.uij[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]) +
#         sum(sm.Lcapapd[j][k][m]*ocuts[o].α10[j,k,m]*mas.ujk[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]) -
#         sum(sm.bigM*mas.uij[i,j,m]*ocuts[o].α12[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]) -
#         sum(sm.bigM*mas.ujk[j,k,m]*ocuts[o].α13[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]) -
#         sum(sm.bigM*mas.ukl[k,l,m]*ocuts[o].α14[k,l,m] for k=1:sm.K for l=1:sm.L for m=1:sm.Mkl[k,l]) <= mas.θ 
#     # dot(objc.*(w*dfo.pi1[o]+(1-w)*dfo.pi2[o]), all_variables(mas.m)[1:end-1]) <= mas.θ            
# end)


# for i=1:5
#     @show file = joinpath(RESULTS_DIR, "sminstances/")*ins[i]
#     fname = file[end-11:end-4]

#     sm = Data(file);
#     lex1 = build_scndmip(1,0,sm.I,sm.J,sm.K,sm.L,sm.Jmax,sm.Kmax,sm.Mij,sm.Mjk,sm.Mkl,sm.demand,sm.bigM,sm.cas,sm.cap,sm.cad,sm.fc,sm.tcp,sm.tcd,sm.tcc,sm.cep,sm.ced,sm.cec,sm.vcs,sm.vcp,sm.vcd,sm.ves,sm.vep,sm.ved,sm.Lcapasp,sm.Lcapapd,sm.gij,sm.gjk,sm.gkl);
#     optimize!(lex1)
#     f1min,f2max = getobjval(lex1,sm)
#     lex2 = build_scndmip(0,1,sm.I,sm.J,sm.K,sm.L,sm.Jmax,sm.Kmax,sm.Mij,sm.Mjk,sm.Mkl,sm.demand,sm.bigM,sm.cas,sm.cap,sm.cad,sm.fc,sm.tcp,sm.tcd,sm.tcc,sm.cep,sm.ced,sm.cec,sm.vcs,sm.vcp,sm.vcd,sm.ves,sm.vep,sm.ved,sm.Lcapasp,sm.Lcapapd,sm.gij,sm.gjk,sm.gkl);
#     optimize!(lex2)
#     f1max,f2min = getobjval(lex2,sm)
#     global fcuts = []; global ocuts = []; global subtime = []; global YN = []; global iter=[];
#     ###Benders Decomposition Algorithm 
#     tCPU = @CPUelapsed solpt,fcut,ocut,subt,iter= DichoBD(sm,f1min,f2max,f1max,f2min);
#     ndp = NDfilter(solpt)

#     tb = DataFrame(fname =[], fcut = [], ocut = [], iter = [], subCPU = [] ,totalCPU=[]);
#     push!(tb.fname, "$fname"); push!(tb.iter, iter);
#     push!(tb.fcut, length(fcut)); push!(tb.ocut, length(ocut)); 
#     push!(tb.subCPU, round(subt; digits=2)); push!(tb.totalCPU, round(tCPU; digits=2));
#     CSV.write(joinpath(RESULTS_DIR, "ClassicBD_tb.csv"), tb, append=true, header=false, delim=' ')

#     otable = zeros(length(ndp),2)
#     for i=1:length(ndp)
#         for j=1:2
#             otable[i,:] = ndp[i]
#         end
#     end
#     CSV.write(joinpath(RESULTS_DIR, "ndp/")*fname*"Y.log",DataFrame(otable, :auto),append=false, header=false, delim=' ' )
# end

