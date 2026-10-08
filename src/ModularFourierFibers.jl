module ModularFourierFibers

using LinearAlgebra, Random
import Nemo
import TOML

export allpairs, observation_set, split_field, sine_field, fourier, pack, unpack,
       predict, target, residual, jacobian, exact_matrix, exact_rank, gauge_matrix,
       kernel_report, block_report, subset_report, rational_memorizer,
       obstruction_report, refine_fiber, search_interpolant, certify_interpolant,
       save_witness, load_witness, write_report, qbilinear, cubic, derivative_rank

const QQ = Nemo.QQ
allpairs(p) = [(a,b) for a in 1:p for b in 1:p]

"Fixed, one-based observation sets; n=22 matches the certified comparison problem."
function observation_set(p, n; seed=1)
    if p == 5 && n in (22,23)
        missing = n == 22 ? [(2,4),(4,4),(5,4)] : [(2,4),(5,4)]
        return setdiff(allpairs(p), missing)
    end
    1 <= n <= p^2 || throw(ArgumentError("need 1 <= n <= p^2"))
    sort(shuffle(MersenneTwister(seed), allpairs(p))[1:n])
end

function split_field(p, start=1_000_000)
    q = Int(Nemo.next_prime(Nemo.ZZ(start)))
    while mod(q-1,4p) != 0; q=Int(Nemo.next_prime(Nemo.ZZ(q))); end
    F=Nemo.GF(q)
    for g in 2:q-1
        z=F(g)^((q-1)÷(4p))
        if z^(2p)!=1 && z^4!=1
            return F,z,q
        end
    end
    error("no primitive 4p-th root found")
end

function sine_field(p)
    R,x=Nemo.polynomial_ring(QQ,"x")
    f=p==5 ? x^4-5x^2+5 : p==7 ? x^6-7x^4+14x^2-7 : error("p=5 or 7 required")
    Nemo.number_field(f,"u")
end

function trig_tables(p; arithmetic=:float, start=1_000_000)
    if arithmetic==:float
        return Float64, cos.(2pi*(0:p-1)/p), sin.(2pi*(0:p-1)/p)
    elseif arithmetic==:finite
        F,z,_=split_field(p,start); w=z^4; ii=z^p
        return F,[(w^a+w^(-a))/2 for a in 0:p-1],[(w^a-w^(-a))/(2ii) for a in 0:p-1]
    elseif arithmetic==:numberfield
        F,u=sine_field(p); t=u^4-4u^2+2; s=p==5 ? u/2 : (u^3-3u)/2
        C=fill(zero(F),p); S=fill(zero(F),p); C[1]=one(F);C[2]=t/2;S[2]=s
        for a in 2:p-1; C[a+1]=t*C[a]-C[a-1];S[a+1]=t*S[a]-S[a-1];end
        return F,C,S
    end
    error("unknown arithmetic")
end

"Explicit exact narrow (2p-1) or wide (4p-3) Fourier realization."
function fourier(p; kind=:narrow, arithmetic=:float, start=1_000_000, phases=nothing)
    isodd(p) && p>=3 || throw(ArgumentError("odd p >= 3 required"))
    kind in (:narrow,:wide) || throw(ArgumentError("unknown realization"))
    F,C,S=trig_tables(p;arithmetic,start); T=eltype(C)
    av=Vector{T}[];bv=Vector{T}[];vv=Vector{T}[];v0=fill(F(1)/p,p)
    for k in 1:(p-1)÷2
        ck=[C[mod(k*a,p)+1] for a in 0:p-1];sk=[S[mod(k*a,p)+1] for a in 0:p-1]
        ca,sa,cb,sb,cv,sv=ck,sk,ck,sk,ck,sk
        if phases!==nothing
            arithmetic==:float || error("phases only supported in floating arithmetic")
            pa,pb=phases[k];ang=2pi*k*(0:p-1)/p
            ca=cos.(ang .+ pa);sa=sin.(ang .+ pa);cb=cos.(ang .+ pb);sb=sin.(ang .+ pb)
            cv=cos.(ang .+ (pa+pb));sv=sin.(ang .+ (pa+pb))
        end
        if kind==:narrow
            for (a,b,v) in [(ca,cb,cv),(sa,-sb,cv),(sa,cb,sv),(ca,sb,sv)]
                push!(av,a);push!(bv,b);push!(vv,v./p)
            end
            v0 .-= 2 .* (cv .+ sv)./p
        else
            for (a,b,v) in [(ca,cb,cv),(sa,sb,-cv),(sa,cb,sv),(ca,sb,sv)]
                push!(av,a);push!(bv,b);push!(vv,v./(2p))
                push!(av,a);push!(bv,-b);push!(vv,-v./(2p))
            end
        end
    end
    push!(av,fill(one(F),p));push!(bv,fill(zero(F),p));push!(vv,v0)
    F,hcat(av...),hcat(bv...),hcat(vv...)
end

pack(A,B,V)=reduce(vcat,[vcat(A[:,j],B[:,j],V[:,j]) for j in axes(A,2)])
function unpack(x,p)
    X=reshape(x,3p,:); X[1:p,:],X[p+1:2p,:],X[2p+1:3p,:]
end
target(p,D,F=Float64)=[F(mod(a+b-2,p)==c-1 ? 1 : 0) for c in 1:p for (a,b) in D]
predict(A,B,V,D)=[sum(V[c,j]*(A[a,j]+B[b,j])^2 for j in axes(A,2)) for c in axes(V,1) for (a,b) in D]
function residual(x,p,D)
    A,B,V=unpack(x,p);predict(A,B,V,D)-target(p,D)
end
function jacobian(A,B,V,D)
    p,K=size(A);n=length(D);J=fill(zero(A[1]),p*n,3p*K)
    for j in 1:K,(i,(a,b)) in enumerate(D),c in 1:p
        r=(c-1)*n+i;o=(j-1)*3p;z=A[a,j]+B[b,j]
        J[r,o+a]=2z*V[c,j];J[r,o+p+b]=2z*V[c,j];J[r,o+2p+c]=z^2
    end
    J
end
exact_matrix(F,A)=Nemo.matrix(F,A)
exact_rank(F,A)=Nemo.rank(exact_matrix(F,A))
function derivative_rank(J;rtol=1e-10)
    S=svd(J;full=true);r=count(>(rtol*S.S[1]),S.S); S,r
end
bigrows(p,D)=[(c-1)*p^2+(a-1)*p+b for c in 1:p for (a,b) in D]

function gauge_matrix(A,B,V,F)
    p,K=size(A);G=Nemo.zero_matrix(F,3p*K,2K)
    for j in 1:K
        o=(j-1)*3p
        for a in 1:p
            G[o+a,2j-1]=1;G[o+p+a,2j-1]=-1
            G[o+a,2j]=A[a,j];G[o+p+a,2j]=B[a,j];G[o+2p+a,2j]=-2V[a,j]
        end
    end
    G
end

function kernel_report(p; arithmetic=:finite,start=1_000_000)
    F,A,B,V=fourier(p;arithmetic,start);J=exact_matrix(F,jacobian(A,B,V,allpairs(p)))
    @assert predict(A,B,V,allpairs(p))==target(p,allpairs(p),F)
    G=gauge_matrix(A,B,V,F);@assert iszero(J*G)
    r=Nemo.rank(J);g=Nemo.rank(G);dims=Int[]
    for js in vcat([collect(4k-3:4k) for k in 1:(p-1)÷2],[[2p-1]])
        cols=reduce(vcat,[collect(3p*(j-1)+1:3p*j) for j in js])
        push!(dims,length(cols)-Nemo.rank(J[:,cols]))
    end
    Dict("p"=>p,"K"=>size(A,2),"arithmetic"=>string(arithmetic),"rank"=>r,
        "parameters"=>Nemo.ncols(J),"kernel"=>Nemo.ncols(J)-r,"gauge"=>g,
        "extra"=>Nemo.ncols(J)-r-g,"localized_kernels"=>dims,
        "coupled_kernel"=>Nemo.ncols(J)-r-sum(dims),"field"=>string(F))
end

function block_report(p;kind=:wide,arithmetic=:finite)
    F,A,B,V=fourier(p;kind,arithmetic);_,C,S=trig_tables(p;arithmetic)
    J=exact_matrix(F,jacobian(A,B,V,allpairs(p)))
    modes=[fill(one(F),p)]
    for k in 1:(p-1)÷2
        push!(modes,[C[mod(k*a,p)+1] for a in 0:p-1]);push!(modes,[S[mod(k*a,p)+1] for a in 0:p-1])
    end
    blocks=[exact_matrix(F,[sum(v[c]*J[(c-1)*p^2+i,j] for c in 1:p) for i in 1:p^2,j in 1:Nemo.ncols(J)]) for v in modes]
    dims=Nemo.rank.(blocks);r=Nemo.rank(J)
    Dict("p"=>p,"K"=>size(A,2),"kind"=>string(kind),"full_rank"=>r,
        "projected_block_dimensions"=>dims,"direct_sum_verified"=>r==sum(dims)),F,J,blocks
end

function subset_report(p;kind=:wide,seeds=1:3,sizes=sort(unique(vcat(collect(5:5:p^2),[p^2]))))
    report,F,J,blocks=block_report(p;kind);rows=Dict[]
    for seed in seeds,n in sizes
        D=observation_set(p,n;seed);rr=Nemo.rank(J[bigrows(p,D),:]);inds=[(a-1)*p+b for (a,b) in D]
        br=[Nemo.rank(U[inds,:]) for U in blocks]
        kind==:wide && @assert rr==sum(br)
        push!(rows,Dict("seed"=>seed,"n"=>n,"rank"=>rr,"naive_rank"=>min(p*n,report["full_rank"]),"block_ranks"=>br))
    end
    report["samples"]=rows;report
end

include("obstructions.jl")
include("numerical.jl")
include("certificate.jl")

function write_report(path, report)
    mkpath(dirname(path));open(path,"w") do io;TOML.print(io,report;sorted=true);end
    println("wrote ",path)
end
end
