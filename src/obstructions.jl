"Polarization of the quadratic Taylor term: Q(h,h)=2 dv*z*dz + v*dz^2."
function qbilinear(A,B,V,D,x,y)
    p,K=size(A);n=length(D);out=fill(zero(A[1]),p*n)
    for j in 1:K,(i,(a,b)) in enumerate(D),c in 1:p
        o=(j-1)*3p;dx=x[o+a]+x[o+p+b];dy=y[o+a]+y[o+p+b];z=A[a,j]+B[b,j]
        out[(c-1)*n+i]+=x[o+2p+c]*z*dy+y[o+2p+c]*z*dx+V[c,j]*dx*dy
    end
    out
end

function cubic(A,D,x,y=x,z=x)
    p,K=size(A);n=length(D);out=fill(zero(A[1]),p*n)
    for j in 1:K,(i,(a,b)) in enumerate(D),c in 1:p
        o=(j-1)*3p;dx=x[o+a]+x[o+p+b];dy=y[o+a]+y[o+p+b];dz=z[o+a]+z[o+p+b]
        out[(c-1)*n+i]+=(x[o+2p+c]*dy*dz+y[o+2p+c]*dx*dz+z[o+2p+c]*dx*dy)/3
    end
    out
end

function kernel_complement(G,N)
    M=G;r=Nemo.rank(G);inds=Int[]
    for j in 1:Nemo.ncols(N)
        next=hcat(M,N[:,j:j]);rr=Nemo.rank(next)
        if rr>r;push!(inds,j);M=next;r=rr;end
    end
    @assert r==Nemo.ncols(N)
    N[:,inds]
end

"Actual coefficients of Q(E*s,E*s), including the factor 2 on off-diagonal monomials."
function quadratic_coefficients(A,B,V,D,F,E)
    k=Nemo.ncols(E);Q=Nemo.zero_matrix(F,size(A,1)*length(D),k*(k+1)÷2);col=0
    for i in 1:k,j in i:k
        col+=1;v=qbilinear(A,B,V,D,E[:,i],E[:,j]);factor=i==j ? 1 : 2
        for row in eachindex(v);Q[row,col]=factor*v[row];end
    end
    Q
end

function cubic_coefficients(A,B,V,D,F,E,Y)
    k=Nemo.ncols(E);pairs=[(i,j) for i in 1:k for j in i:k];idx=Dict(t=>i for (i,t) in enumerate(pairs))
    triples=[(i,j,l) for i in 1:k for j in i:k for l in j:k]
    T=Nemo.zero_matrix(F,size(A,1)*length(D),length(triples))
    for (col,(i,j,l)) in enumerate(triples)
        ei=E[:,i];ej=E[:,j];el=E[:,l]
        if i==l
            corr=2 .* qbilinear(A,B,V,D,ei,Y[:,idx[(i,i)]])
            raw=cubic(A,D,ei)
        elseif i==j
            corr=2 .* (qbilinear(A,B,V,D,el,Y[:,idx[(i,i)]]) .+ qbilinear(A,B,V,D,ei,Y[:,idx[(i,l)]]))
            raw=3 .* cubic(A,D,ei,ei,el)
        elseif j==l
            corr=2 .* (qbilinear(A,B,V,D,ei,Y[:,idx[(j,j)]]) .+ qbilinear(A,B,V,D,ej,Y[:,idx[(i,j)]]))
            raw=3 .* cubic(A,D,ei,ej,ej)
        else
            corr=2 .* (qbilinear(A,B,V,D,ei,Y[:,idx[(j,l)]]) .+ qbilinear(A,B,V,D,ej,Y[:,idx[(i,l)]]) .+ qbilinear(A,B,V,D,el,Y[:,idx[(i,j)]]))
            raw=6 .* cubic(A,D,ei,ej,el)
        end
        for row in 1:Nemo.nrows(T);T[row,col]=corr[row]+raw[row];end
    end
    T
end

function obstruction_report(n; arithmetic=:finite,start=1_000_000)
    p=5;F,A,B,V=fourier(p;arithmetic,start);D=observation_set(p,n);all=allpairs(p)
    JD=exact_matrix(F,jacobian(A,B,V,D));JA=exact_matrix(F,jacobian(A,B,V,all))
    kd,N=Nemo.nullspace(JD);@assert Nemo.rank(JD)==Nemo.rank(JA) && iszero(JA*N)
    G=gauge_matrix(A,B,V,F);E=kernel_complement(G,N)
    _,LD=Nemo.nullspace(transpose(JD));_,LA=Nemo.nullspace(transpose(JA))
    QD=quadratic_coefficients(A,B,V,D,F,E);QA=quadratic_coefficients(A,B,V,all,F,E)
    @assert iszero(transpose(LD)*QD) && iszero(transpose(LA)*QA)
    ok,Y=Nemo.can_solve_with_solution(JD,-QD;side=:right);@assert ok && iszero(JA*Y+QA)
    TD=cubic_coefficients(A,B,V,D,F,E,Y);TA=cubic_coefficients(A,B,V,all,F,E,Y)
    rd=Nemo.rank(transpose(LD)*TD);ra=Nemo.rank(transpose(LA)*TA)
    Dict("p"=>p,"K"=>size(A,2),"n"=>n,"arithmetic"=>string(arithmetic),"field"=>string(F),
        "training_rank"=>Nemo.rank(JD),"full_rank"=>Nemo.rank(JA),"kernel"=>kd,
        "gauge"=>Nemo.rank(G),"slice_dimension"=>Nemo.ncols(E),"quadratic_monomials"=>Nemo.ncols(QD),
        "quadratic_obstruction_rank_D"=>0,"quadratic_obstruction_rank_all"=>0,
        "cubic_monomials"=>Nemo.ncols(TD),"cubic_obstruction_rank_D"=>rd,"cubic_obstruction_rank_all"=>ra)
end
