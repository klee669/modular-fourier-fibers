using ModularFourierFibers, LinearAlgebra, Random
import Nemo

BLAS.set_num_threads(1)
const ROOT=dirname(@__DIR__)
const OUT=get(ENV,"MFF_RESULTS",joinpath(ROOT,"results"))
emit(name,r)=write_report(joinpath(OUT,name*".toml"),r)

"One presentation report: independent output-frequency blocks and their restrictions."
function eight_neuron_result()
    println("[1/3] Eight neurons per frequency pair: output-frequency blocks")
    rows=Dict[]
    for p in (5,7)
        r,F,_,_=block_report(p;kind=:eight_per_frequency,arithmetic=:numberfield)
        @assert r["direct_sum_verified"]
        @assert r["projected_block_dimensions"]==vcat([4p-3],fill(8p-15,p-1))
        r["field"]=string(F);push!(rows,r)
        println("  p=$p, K=$(r["K"]): rank=$(r["full_rank"]), blocks=$(r["projected_block_dimensions"])")
    end
    sweep=subset_report(7;kind=:eight_per_frequency,arithmetic=:numberfield,seeds=1:1,sizes=[20,25,30,35,40,45,49])
    emit("eight_neurons_per_frequency",Dict(
        "construction"=>"8 neurons per conjugate frequency pair plus one constant neuron; K=4p-3",
        "evidence"=>"exact characteristic-zero cyclotomic number-field linear algebra for the listed cases",
        "blocks"=>rows,"partial_observations_p7"=>sweep,
        "limitation"=>"No general-p theorem or nonlinear local completion claim for this construction."))
end

"One presentation report: kernel decomposition and common partial/full tangent kernels."
function four_neuron_result()
    println("[2/3] Four neurons per frequency pair: kernels and partial observations")
    rows=Dict[]
    for p in (5,7)
        r=kernel_report(p;arithmetic=:numberfield)
        @assert r["rank"]==(p==5 ? 105 : 229)
        @assert r["gauge"]==4p-2 && r["extra"]==3p-3
        push!(rows,r)
        println("  p=$p, K=$(r["K"]): kernel=$(r["kernel"]) = gauge $(r["gauge"]) + extra $(r["extra"])")
    end
    p=5;F,A,B,V=fourier(p;kind=:four_per_frequency,arithmetic=:numberfield)
    JA=exact_matrix(F,jacobian(A,B,V,allpairs(p)));ra=Nemo.rank(JA);comparisons=Dict[]
    for n in (22,23)
        D=observation_set(p,n);JD=exact_matrix(F,jacobian(A,B,V,D));rd=Nemo.rank(JD)
        # The full system contains all partial equations. Equal ranks then imply equal kernels.
        @assert rd==ra==105
        push!(comparisons,Dict("p"=>p,"K"=>9,"n"=>n,"D"=>[[a,b] for (a,b) in D],
            "arithmetic"=>"numberfield","field"=>string(F),"training_rank"=>rd,"full_rank"=>ra,
            "common_kernel_dimension"=>135-rd,"kernels_equal"=>true))
        println("  |D|=$n: training rank = full rank = $ra; tangent kernels coincide")
    end
    emit("four_neurons_per_frequency",Dict(
        "construction"=>"4 neurons per conjugate frequency pair plus one constant neuron; K=2p-1",
        "evidence"=>"exact characteristic-zero number-field linear algebra",
        "kernels"=>rows,"partial_full_comparisons"=>comparisons,
        "limitation"=>"Tangent-kernel equality does not prove smoothness or equality of fiber germs."))
end

function ranks(;full=false)
    for p in (5,7)
        r=kernel_report(p);@assert r["rank"]==(p==5 ? 105 : 229)
        @assert r["gauge"]==4p-2 && r["extra"]==3p-3
        emit("kernel_p$p",r)
        emit("kernel_p$(p)_second_prime",kernel_report(p;start=2_000_000))
        full && emit("kernel_p$(p)_numberfield",kernel_report(p;arithmetic=:numberfield))
    end
    for p in (full ? (3,5,7) : (5,7))
        r,_,_,_=block_report(p);@assert r["direct_sum_verified"]
        emit("eight_per_frequency_blocks_p$p",r)
        nr,_,_,_=block_report(p;kind=:four_per_frequency);emit("four_per_frequency_projected_blocks_p$p",nr)
    end
    for p in (5,7)
        emit("eight_per_frequency_sweep_p$p",subset_report(p))
        emit("four_per_frequency_sweep_p$p",subset_report(p;kind=:four_per_frequency,sizes=p==5 ? [19,20,21,22,23,24,25] : [30,33,34,35,36,40,49]))
    end
end

function moderate()
    for p in (5,7)
        K=4p-3;n=round(Int,0.4p^2);D=observation_set(p,n)
        A,B,V=rational_memorizer(p,K,D)
        F,AF,BF,VF=fourier(p;kind=:eight_per_frequency,arithmetic=:finite)
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
    println("[3/3] Nearby completion: Newton correction using only observed pairs")
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
    trainmax=maximum(r["training_residual"] for r in rows)
    fullmax=maximum(r["full_residual"] for r in rows)
    println("  $(length(rows)) samples: max training residual=$trainmax; max full residual=$fullmax")
    emit("nearby_completion",Dict("p"=>5,"K"=>9,"construction"=>"4 neurons per frequency pair",
        "sample_count"=>length(rows),"seed"=>20260927,
        "max_training_residual"=>trainmax,"max_full_residual"=>fullmax,
        "samples"=>rows,"status"=>"numerical evidence only; not a germ-equality certificate"))
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

command=isempty(ARGS) ? "showcase" : ARGS[1]
if command in ("showcase","quick","full")
    # Default and compatibility aliases deliberately produce only these three reports.
    eight_neuron_result();four_neuron_result();nearby()
elseif command=="eight";eight_neuron_result()
elseif command=="four";four_neuron_result()
elseif command in ("ranks","rank");ranks(;full=true)
elseif command=="jets";jets()
elseif command=="moderate";moderate()
elseif command=="realizations";realizations()
elseif command=="nearby";nearby()
elseif command=="search";discover()
elseif command=="certify";certify()
else
    error("usage: julia --project=. scripts/run.jl [showcase|full|quick|eight|four|nearby|ranks|jets|moderate|realizations|search|certify [witness.txt]]")
end
