function smSCND_MIP(sm,tolerance)#sm.I,J,K,L,Jmax,Kmax,Mij,Mjk,Mkl,demand,bigM,,sm.ves,sm.vep,sm.ved,sm.Lcapasp,sm.Lcapapd,sm.gij,sm.gjk,sm.gkl
    scnd = vModel(CPLEX.Optimizer);
    MOI.set(scnd, MOI.RelativeGapTolerance(), tolerance)
    set_silent(scnd)
    @variable(scnd, y[1:sm.J+sm.K,1:2], Bin)
    @variable(scnd, uij[i=1:sm.I,j=1:sm.J,1:sm.Mij[i,j]], Bin)
    @variable(scnd, ujk[j=1:sm.J,k=1:sm.K,1:sm.Mjk[j,k]], Bin)
    @variable(scnd, ukl[k=1:sm.K,l=1:sm.L,1:sm.Mkl[k,l]], Bin)
    @variable(scnd, 0<= xij[i=1:sm.I,j=1:sm.J,1:sm.Mij[i,j]])
    @variable(scnd, 0<= xjk[j=1:sm.J,k=1:sm.K,1:sm.Mjk[j,k]])
    @variable(scnd, 0<= xkl[k=1:sm.K,l=1:sm.L,1:sm.Mkl[k,l]])
    @variable(scnd, 0<= h[1:sm.J+sm.K,1:2])
    
    @addobjective(scnd, Min, sum(sm.fc[j,t]*y[j,t] for j=1:sm.J+sm.K for t=1:2) +
        sum(sm.gij[i][j][m]*uij[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j])+
        sum(sm.gjk[j][k][m]*ujk[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]) +
        sum(sm.gkl[k][l][m]*ukl[k,l,m] for k=1:sm.K for l=1:sm.L for m=1:sm.Mkl[k,l]) +
        sum(sm.vcs[i]*xij[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]) +
        sum([sm.vcp;sm.vcd][j,t]*h[j,t] for j=1:sm.J+sm.K for t=1:2) +
        sum(sm.tcp[i][j][m]*xij[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]) +
        sum(sm.tcd[j][k][m]*xjk[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]) +
        sum(sm.tcc[k][l][m]*xkl[k,l,m] for k=1:sm.K for l=1:sm.L for m=1:sm.Mkl[k,l]) )
    @addobjective(scnd, Min, sum(sm.ves[i]*xij[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]) +
        sum([sm.vep;sm.ved][j,t]*h[j,t] for j=1:sm.J+sm.K for t=1:2) +
        sum(sm.cep[i][j][m]*xij[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j] )+
        sum(sm.ced[j][k][m]*xjk[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k] )+
        sum(sm.cec[k][l][m]*xkl[k,l,m] for k=1:sm.K for l=1:sm.L for m=1:sm.Mkl[k,l] ) );
    ######### constraint 3 #############
    @constraints(scnd, begin
        [j=1:sm.J], sum(xij[i,j,m] for i=1:sm.I for m=1:sm.Mij[i,j]) == sum(xjk[j,k,m] for k=1:sm.K for m=1:sm.Mjk[j,k]);
        [k=1:sm.K], sum(xjk[j,k,m] for j=1:sm.J for m=1:sm.Mjk[j,k]) == sum(xkl[k,l,m] for l=1:sm.L for m=1:sm.Mkl[k,l]);
    end)
    ########### constraint 4-6 #############
    @constraints(scnd, begin
        [j=1:sm.J], sum(h[j,:]) == sum(xij[i,j,m] for i=1:sm.I for m=1:sm.Mij[i,j]);
        [k=1:sm.K], sum(h[k+sm.J,:] ) == sum(xjk[j,k,m] for j=1:sm.J for m=1:sm.Mjk[j,k]);
        [l=1:sm.L], sum(xkl[k,l,m] for k=1:sm.K for m=1:sm.Mkl[k,l]) >= sm.demand[l];
    end );
    ########### constraint 7-9 #############
    @constraint(scnd,con7[i=1:sm.I], sum(xij[i,j,m] for j=1:sm.J for m=1:sm.Mij[i,j] ) <= sm.cas[i]);
    @constraint(scnd,con81[j=1:sm.J, t=1:2], sum(h[j,t]) <= sm.cap[j]*y[j,t]);
    @constraint(scnd,con82[k=1:sm.K, t=1:2], sum(h[sm.J+k,t]) <= sm.cad[k]*y[sm.J+k,t]);
    @constraint(scnd,con9[j=1:sm.J+sm.K], sum(y[j,:]) <= 1)
    # ########### constraint 10 ############# 
    # @constraint(scnd,con101[i=1:sm.I,j=1:sm.J], sum(uij[i,j,m] for m=1:sm.Mij[i,j]) <= 1); #Those constraints from the paper provide incorrect solutions
    # @constraint(scnd,con102[j=1:sm.J,k=1:sm.K], sum(ujk[j,k,m] for m=1:sm.Mjk[j,k]) <= 1);
    # @constraint(scnd,con103[k=1:sm.K,l=1:sm.L], sum(ukl[k,l,m] for m=1:sm.Mkl[k,l]) <= 1);
    @constraint(scnd,[i=1:sm.I,j=1:sm.J], sum(uij[i,j,m] for m=1:sm.Mij[i,j]) <=  sum(y[j,:]));
    @constraint(scnd,[j=1:sm.J,k=1:sm.K], sum(ujk[j,k,m] for m=1:sm.Mjk[j,k]) <= (sum(y[j,:])+sum(y[sm.J+k,:]))/2 );
    @constraint(scnd,[k=1:sm.K,l=1:sm.L], sum(ukl[k,l,m] for m=1:sm.Mkl[k,l]) <= sum(y[sm.J+k,:]));
    # ########### constraint 11 ############# This includes redundant constraints: where there is just one mode available between nodes, i.e. Mij[i,j] = 1.
    # @constraint(scnd,con121[i=1:sm.I, j=1:sm.J, m=1:sm.Mij[i,j]], xij[i,j,m] >= sm.Lcapasp[i][j][m]*uij[i,j,m] );
    # @constraint(scnd,con122[j=1:sm.J, k=1:sm.K, m=1:sm.Mjk[j,k]], xjk[j,k,m] >= sm.Lcapapd[j][k][m]*ujk[j,k,m]);
    # Thus, the below is the new version
    for i=1:sm.I
        for j=1:sm.J
            if sm.Mij[i,j] == 2
                @constraint(scnd, xij[i,j,2] >= sm.Lcapasp[i][j][2]*uij[i,j,2] );
            end
        end
    end            
    for j=1:sm.J
        for k=1:sm.K
            if sm.Mjk[j,k] == 2
                @constraint(scnd, xjk[j,k,2] >= sm.Lcapapd[j][k][2]*ujk[j,k,2]);
            end
        end
    end  
    ########### constraint 12 #############
    @constraint(scnd, con12[i=1:sm.I,j=1:sm.J,m=1:sm.Mij[i,j]], xij[i,j,m] <= sm.bigM*uij[i,j,m]);
    @constraint(scnd, con13[j=1:sm.J,k=1:sm.K,m=1:sm.Mjk[j,k]], xjk[j,k,m] <= sm.bigM*ujk[j,k,m]);
    @constraint(scnd, con14[k=1:sm.K, l=1:sm.L, m=1:sm.Mkl[k,l]], xkl[k,l,m] <= sm.bigM*ukl[k,l,m]);
    ########### constraint 13-14 #############
    @constraint(scnd, sum(y[j,t] for j=1:sm.J for t=1:2) <= sm.Jmax);
    @constraint(scnd, sum(y[k,t] for k=sm.J+1:sm.K+sm.J for t=1:2) <= sm.Kmax)
    return scnd
end

# function getobjval(model,sm)
#     y = value.(model[:y])
#     uij = value.(model[:uij])
#     ujk = value.(model[:ujk])
#     ukl = value.(model[:ukl])
#     xij = value.(model[:xij])
#     xjk = value.(model[:xjk])
#     xkl = value.(model[:xkl])
#     h = value.(model[:h])

#     obj1 = sum(sm.fc[j,t]*y[j,t] for j=1:sm.J+sm.K for t=1:2) +
#         sum(sm.gij[i][j][m]*uij[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j])+
#         sum(sm.gjk[j][k][m]*ujk[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]) +
#         sum(sm.gkl[k][l][m]*ukl[k,l,m] for k=1:sm.K for l=1:sm.L for m=1:sm.Mkl[k,l]) +
#         sum(sm.vcs[i]*xij[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]) +
#         sum([sm.vcp;sm.vcd][j,t]*h[j,t] for j=1:sm.J+sm.K for t=1:2) +
#         sum(sm.tcp[i][j][m]*xij[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]) +
#         sum(sm.tcd[j][k][m]*xjk[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]) +
#         sum(sm.tcc[k][l][m]*xkl[k,l,m] for k=1:sm.K for l=1:sm.L for m=1:sm.Mkl[k,l]) 
        
#     obj2 = (sum(sm.ves[i]*xij[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j]) +
#     sum([sm.vep;sm.ved][j,t]*h[j,t] for j=1:sm.J+sm.K for t=1:2) +
#     sum(sm.cep[i][j][m]*xij[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j] )+
#     sum(sm.ced[j][k][m]*xjk[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k] )+
#     sum(sm.cec[k][l][m]*xkl[k,l,m] for k=1:sm.K for l=1:sm.L for m=1:sm.Mkl[k,l] ))
#     # obj1 = sum(sm.fc[j,t]*y[j,t] for j=1:sm.J+sm.K for t=1:2) +
#     # sum(sm.gij[i][j][m]*uij[i,j,m] for i=1:sm.I for j=1:sm.J for m=1:sm.Mij[i,j])+
#     # sum(sm.gjk[j][k][m]*ujk[j,k,m] for j=1:sm.J for k=1:sm.K for m=1:sm.Mjk[j,k]) +
#     # sum(sm.gkl[k][l][m]*ukl[k,l,m] for k=1:sm.K for l=1:sm.L for m=1:sm.Mkl[k,l]) +
#     # sum(sm.vcs[i]*xij[i,j,1] for i=1:sm.I for j=1:sm.J ) + sum(sm.vcs[i]*xij[i,j,2] for i=1:sm.I for j=1:sm.J if sm.Mij[i,j]==2) +
#     # sum([sm.vcp;sm.vcd][j,t]*h[j,t] for j=1:sm.J+sm.K for t=1:2) +
#     # sum(sm.tcp[i][j][1]*xij[i,j,1] for i=1:sm.I for j=1:sm.J) + sum(sm.tcp[i][j][2]*xij[i,j,2] for i=1:sm.I for j=1:sm.J if sm.Mij[i,j]==2) +
#     # sum(sm.tcd[j][k][1]*xjk[j,k,1] for j=1:sm.J for k=1:sm.K) + sum(sm.tcd[j][k][2]*xjk[j,k,2] for j=1:sm.J for k=1:sm.K if sm.Mjk[j,k]==2) +
#     # sum(sm.tcc[k][l][m]*xkl[k,l,m] for k=1:sm.K for l=1:sm.L for m=1:sm.Mkl[k,l]) 
    
#     # obj2 = sum(sm.ves[i]*xij[i,j,1] for i=1:sm.I for j=1:sm.J) + sum(sm.ves[i]*xij[i,j,2] for i=1:sm.I for j=1:sm.J if sm.Mij[i,j]==2) +
#     # sum([sm.vep;sm.ved][j,t]*h[j,t] for j=1:sm.J+sm.K for t=1:2) +
#     # sum(sm.cep[i][j][1]*xij[i,j,1] for i=1:sm.I for j=1:sm.J)+ sum(sm.cep[i][j][2]*xij[i,j,2] for i=1:sm.I for j=1:sm.J if sm.Mij[i,j]==2 )+
#     # sum(sm.ced[j][k][1]*xjk[j,k,1] for j=1:sm.J for k=1:sm.K)+ sum(sm.ced[j][k][2]*xjk[j,k,2] for j=1:sm.J for k=1:sm.K if sm.Mjk[j,k]==2 )+
#     # sum(sm.cec[k][l][m]*xkl[k,l,m] for k=1:sm.K for l=1:sm.L for m=1:sm.Mkl[k,l] )
#     return round(obj1),round(obj2)#obj1,obj2       
# end
