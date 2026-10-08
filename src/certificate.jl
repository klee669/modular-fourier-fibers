to_qq(x::AbstractFloat)=QQ(Rational{BigInt}(x))
to_qq(x)=QQ(x)

"Portable text artifact: exact rationals, one parameter per line, neuron-major alpha/beta/v."
function save_witness(path,x)
    mkpath(dirname(path))
    open(path,"w") do io
        println(io,"# p=5 K=9; neuron-major alpha[1:5], beta[1:5], v[1:5]")
        for xx in x
            q=to_qq(xx);println(io,Nemo.numerator(q),"/",Nemo.denominator(q))
        end
    end
end
function load_witness(path)
    [begin
        n,d=split(strip(line),'/');QQ(parse(BigInt,n))/QQ(parse(BigInt,d))
     end for line in readlines(path) if !isempty(strip(line)) && !startswith(strip(line),"#")]
end

"Exact-rational contraction certificate; floating QR chooses columns but proves nothing."
function certify_interpolant(x; columns=nothing)
    p=5;K=9;D=observation_set(p,22);length(x)==3p*K || error("wrong witness size")
    xr=to_qq.(x);A,B,V=unpack(xr,p);J=exact_matrix(QQ,jacobian(A,B,V,D));m=p*length(D)
    cols=columns===nothing ? Vector(qr(Float64.(jacobian(A,B,V,D)),ColumnNorm()).p)[1:m] : columns
    length(cols)==m && length(unique(cols))==m || error("invalid columns")
    JS=J[:,cols];BI=inv(JS);@assert isone(BI*JS)
    rr=predict(A,B,V,D)-target(p,D,QQ);r=Nemo.matrix(QQ,m,1,rr)
    eta=maximum(abs,(BI*r)[:,1]);eta>0 || error("witness is already exact; choose a nonzero certificate radius")
    bn=maximum(sum(abs(BI[i,j]) for j in 1:m) for i in 1:m);rho=2eta;M=maximum(abs,xr)+rho
    L=24K*M;q=bn*L*rho
    @assert q<1 && eta+q*rho<rho "contraction test failed; refine numerical candidate"
    ra=predict(A,B,V,allpairs(p))-target(p,allpairs(p),QQ)
    omitted=findall(i->!(allpairs(p)[mod(i-1,p^2)+1] in D),eachindex(ra))
    idx=omitted[argmax(abs.(ra[omitted]))];lower=abs(ra[idx])-12K*M^2*rho
    @assert lower>0 "nongeneralization not certified"
    pair=allpairs(p)[mod(idx-1,p^2)+1];c=(idx-1)÷p^2+1
    Dict("p"=>p,"K"=>K,"n"=>length(D),"selected_columns"=>cols,"root_exists"=>true,
        "rank_at_root"=>m,"fiber_dimension"=>length(x)-m,"heldout_coordinate"=>[pair[1],pair[2],c],
        "radius_approx"=>Float64(rho),"contraction_approx"=>Float64(q),"heldout_error_lower_approx"=>Float64(lower),
        "eta_exact"=>string(eta),"radius_exact"=>string(rho),"inverse_inf_norm_exact"=>string(bn),
        "jacobian_lipschitz_bound_exact"=>string(L),"contraction_exact"=>string(q),"heldout_error_lower_exact"=>string(lower),
        "loss_convention"=>"uncentered sum of squared residuals; positive smooth local prior",
        "local_learning_coefficient"=>m/2)
end
