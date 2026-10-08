function rational_memorizer(p,K,D;seed=101)
    rng=MersenneTwister(seed);A=QQ.(rand(rng,-9:9,p,K));B=QQ.(rand(rng,-9:9,p,K));V=QQ.(rand(rng,-9:9,p,K))
    H=exact_matrix(QQ,[(A[a,j]+B[b,j])^2 for (a,b) in D,j in 1:K])
    r,R=Nemo.rref(H);r==length(D) || error("features not full row rank; use another seed or larger K")
    cols=[findfirst(j->!iszero(R[i,j]),1:K) for i in 1:r]
    Y=exact_matrix(QQ,[QQ(mod(a+b-2,p)==c-1 ? 1 : 0) for (a,b) in D,c in 1:p])
    delta=Nemo.solve(H[:,cols],Y-H*exact_matrix(QQ,permutedims(V));side=:right)
    for i in 1:r,c in 1:p;V[c,cols[i]]+=delta[i,c];end
    @assert predict(A,B,V,D)==target(p,D,QQ)
    A,B,V
end

function refine_fiber(x,p,D;maxit=40,atol=2e-13)
    for it in 1:maxit
        r=residual(x,p,D);norm(r,Inf)<atol && return x,it
        A,B,V=unpack(x,p);S,k=derivative_rank(jacobian(A,B,V,D);rtol=1e-12)
        delta=-S.V[:,1:k]*((S.U[:,1:k]'*r)./S.S[1:k]);step=1.0
        while norm(residual(x+step*delta,p,D))>norm(r) && step>1e-5;step/=2;end
        x=x+step*delta
    end
    x,maxit
end

function lmfit(x,p,D;maxit=2000)
    mu=1e-3;r=residual(x,p,D);val=sum(abs2,r)
    for it in 1:maxit
        norm(r,Inf)<1e-12 && return x,it
        A,B,V=unpack(x,p);J=jacobian(A,B,V,D);JJ=J*J';accepted=false
        for trial in 1:15
            dx=-J'*((JJ+mu*I)\r);xn=x+dx;rn=residual(xn,p,D);vn=sum(abs2,rn)
            if vn<val
                x=xn;r=rn;val=vn;mu=max(mu/3,1e-16);accepted=true;break
            end
            mu*=10
        end
        !accepted && return x,it
    end
    x,maxit
end

"Numerical discovery only. Call certify_interpolant to prove existence."
function search_interpolant(;seed=1729,attempts=60)
    p=5;D=observation_set(p,22);_,A,B,V=fourier(p);x0=pack(A,B,V);rng=MersenneTwister(seed)
    for attempt in 1:attempts
        x=attempt<=30 ? 0.5randn(rng,length(x0)) : x0+(0.2+0.1*(attempt-31))*randn(rng,length(x0))
        sol,it=lmfit(x,p,D);AA,BB,VV=unpack(sol,p);s=svdvals(jacobian(AA,BB,VV,D))
        rd=norm(residual(sol,p,D),Inf);ra=norm(residual(sol,p,allpairs(p)),Inf)
        println("attempt=$attempt iterations=$it train=$rd full=$ra smallest_sv=$(s[end])");flush(stdout)
        if rd<1e-10 && ra>1e-5 && s[end]>1e-6
            return sol,Dict("seed"=>seed,"attempt"=>attempt,"iterations"=>it,"train_residual"=>rd,
                "full_residual"=>ra,"smallest_singular_value"=>s[end])
        end
    end
    error("no candidate found in $attempts attempts; no nonexistence conclusion follows")
end
