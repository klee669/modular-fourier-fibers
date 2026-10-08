using ModularFourierFibers, LinearAlgebra, Random
import Nemo

BLAS.set_num_threads(1)
const ROOT=dirname(@__DIR__)
const OUT=get(ENV,"MFF_RESULTS",joinpath(ROOT,"results"))
emit(name,r)=write_report(joinpath(OUT,name*".toml"),r)

function ranks(;full=false)
    for p in (5,7)
        r=kernel_report(p);@assert r["rank"]==(p==5 ? 105 : 229)
        @assert r["gauge"]==4p-2 && r["extra"]==3p-3
        emit("kernel_p$p",r)
        emit("kernel_p$(p)_second_prime",kernel_report(p;start=2_000_000))
        full && emit("kernel_p$(p)_numberfield",kernel_report(p;arithmetic=:numberfield))
    end
    for p in (full ? (3,5,7,11) : (5,7))
        r,_,_,_=block_report(p);@assert r["direct_sum_verified"]
        emit("wide_blocks_p$p",r)
        nr,_,_,_=block_report(p;kind=:narrow);emit("narrow_projected_blocks_p$p",nr)
    end
    for p in (5,7)
        emit("wide_sweep_p$p",subset_report(p))
        emit("narrow_sweep_p$p",subset_report(p;kind=:narrow,sizes=p==5 ? [19,20,21,22,23,24,25] : [30,33,34,35,36,40,49]))
    end
end

function moderate()
    for p in (5,7)
        K=4p-3;n=round(Int,0.4p^2);D=observation_set(p,n)
        A,B,V=rational_memorizer(p,K,D)
        F,AF,BF,VF=fourier(p;kind=:wide,arithmetic=:finite)
        rf=exact_rank(F,jacobian(AF,BF,VF,D));rfa=exact_rank(F,jacobian(AF,BF,VF,allpairs(p)))
        rm=exact_rank(F,jacobian(F.(A),F.(B),F.(V),D));rma=exact_rank(F,jacobian(F.(A),F.(B),F.(V),allpairs(p)))
        @assert rf==rm==p*n && rma==p^3
        @assert predict(A,B,V,allpairs(p))!=target(p,allpairs(p),Nemo.QQ)
        emit("moderate_p$p",Dict("p"=>p,"K"=>K,"n"=>n,"D"=>[[a,b] for (a,b) in D],
            "Fourier_training_rank"=>rf,"memorizer_training_rank"=>rm,
            "Fourier_full_rank"=>rfa,"memorizer_full_rank"=>rma,"LLC_each"=>p*n/2))
    end
end

function jets()
    for n in (22,23),start in (1_000_000,2_000_000)
        r=obstruction_report(n;start)
        @assert r["cubic_obstruction_rank_D"]==r["cubic_obstruction_rank_all"]==0
        emit("jets_n$(n)_$(start)",r)
    end
end

function realizations()
    p=5;D=observation_set(p,22);rows=Dict[]
    for phases in [[(0.13,0.27),(-0.41,0.19)],[(0.7,-0.2),(1.1,0.6)],[(1.4,0.9),(-0.8,-1.2)]]
        _,A,B,V=fourier(p;phases);x=pack(A,B,V);err=norm(residual(x,p,allpairs(p)),Inf)
        _,r=derivative_rank(jacobian(A,B,V,allpairs(p)));_,rd=derivative_rank(jacobian(A,B,V,D))
        @assert err<1e-12 && r==rd==105
        push!(rows,Dict("phases"=>[[a,b] for (a,b) in phases],"full_residual"=>err,"full_rank"=>r,"training_rank"=>rd))
    end
    F,A,B,V=fourier(p;arithmetic=:finite)
    A=hcat(A,A[:,1]);B=hcat(B,B[:,1]);V=hcat(V,V[:,1]/2);V[:,1]./=2
    @assert predict(A,B,V,allpairs(p))==target(p,allpairs(p),F)
    r=exact_rank(F,jacobian(A,B,V,allpairs(p)));@assert r==105
    emit("realizations",Dict("phase_variants"=>rows,"split_K"=>10,"split_rank"=>r,"split_kernel"=>150-r))
end

function nearby()
    p=5;_,A,B,V=fourier(p);x0=pack(A,B,V);rng=MersenneTwister(20260927);rows=Dict[]
    for n in (22,23)
        D=observation_set(p,n);S,r=derivative_rank(jacobian(A,B,V,D));N=S.V[:,r+1:end]
        for scale in (1e-2,1e-3,1e-4),rep in 1:8
            u=N*randn(rng,size(N,2));u/=norm(u)
            sol,it=refine_fiber(x0+scale*u,p,D)
            push!(rows,Dict("n"=>n,"scale"=>scale,"rep"=>rep,"distance"=>norm(sol-x0),
                "training_residual"=>norm(residual(sol,p,D),Inf),"full_residual"=>norm(residual(sol,p,allpairs(p)),Inf),"iterations"=>it))
        end
    end
    emit("nearby_numerical",Dict("samples"=>rows,"status"=>"numerical evidence only; not a germ-equality certificate"))
end

function discover()
    x,r=search_interpolant();cert=certify_interpolant(x)
    save_witness(joinpath(OUT,"discovered_witness.txt"),x)
    emit("search",r);emit("discovered_certificate",cert)
end

function certify()
    path=length(ARGS)>=2 ? abspath(ARGS[2]) : joinpath(ROOT,"data","nongeneralizing_witness.txt")
    r=certify_interpolant(load_witness(path));emit("certificate",r)
    println("CERTIFIED: real interpolant exists; rank=",r["rank_at_root"],
        "; held-out error > ",r["heldout_error_lower_approx"],"; LLC=",r["local_learning_coefficient"])
end

command=isempty(ARGS) ? "quick" : ARGS[1]
if command in ("quick","full")
    ranks(;full=command=="full");moderate();jets();realizations();nearby();certify()
elseif command=="ranks";ranks(;full=true)
elseif command=="jets";jets()
elseif command=="moderate";moderate()
elseif command=="realizations";realizations()
elseif command=="nearby";nearby()
elseif command=="search";discover()
elseif command=="certify";certify()
else
    error("usage: julia --project=. scripts/run.jl [quick|full|ranks|jets|moderate|realizations|nearby|search|certify [witness.txt]]")
end
