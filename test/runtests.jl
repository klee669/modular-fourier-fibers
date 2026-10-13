using Test, ModularFourierFibers, LinearAlgebra, Random
import Nemo
const MF=ModularFourierFibers
BLAS.set_num_threads(1)

@testset "Fourier roots and gauges" begin
    for p in (5,7),kind in (:four_per_frequency,:eight_per_frequency)
        F,A,B,V=fourier(p;kind,arithmetic=:finite);D=allpairs(p)
        @test predict(A,B,V,D)==target(p,D,F)
        G=gauge_matrix(A,B,V,F);J=exact_matrix(F,jacobian(A,B,V,D))
        @test iszero(J*G)
        @test Nemo.rank(G)==2size(A,2)
    end
end

@testset "Taylor identity and polarization (regression for missing base z factor)" begin
    F,A,B,V=fourier(5;arithmetic=:finite);rng=MersenneTwister(44);p=5
    # A random base point, not only a Fourier root, catches cancellation-specific errors.
    A=F.(rand(rng,-4:4,size(A)));B=F.(rand(rng,-4:4,size(B)));V=F.(rand(rng,-4:4,size(V)))
    x=pack(A,B,V);h=F.(rand(rng,-3:3,length(x)));g=F.(rand(rng,-3:3,length(x)))
    D=observation_set(p,22);t=F(7);AA,BB,VV=unpack(x+t.*h,p)
    lhs=predict(AA,BB,VV,D)-predict(A,B,V,D)
    J=jacobian(A,B,V,D)
    @test lhs==t.*(J*h)+t^2 .* qbilinear(A,B,V,D,h,h)+t^3 .* cubic(A,D,h)
    @test qbilinear(A,B,V,D,h+g,h+g)==qbilinear(A,B,V,D,h,h)+2qbilinear(A,B,V,D,h,g)+qbilinear(A,B,V,D,g,g)
    # Check the complete monomial-table assembly against direct substitution.
    k=3;E=exact_matrix(F,F.(rand(rng,-3:3,length(x),k)));Y=exact_matrix(F,F.(rand(rng,-3:3,length(x),6)))
    s=F.([2,3,5]);m2=[s[i]*s[j] for i in 1:k for j in i:k];m3=[s[i]*s[j]*s[l] for i in 1:k for j in i:k for l in j:k]
    es=E*Nemo.matrix(F,k,1,s);ym=Y*Nemo.matrix(F,6,1,m2)
    Q=MF.quadratic_coefficients(A,B,V,D,F,E)
    @test Q*Nemo.matrix(F,6,1,m2)==Nemo.matrix(F,p*length(D),1,qbilinear(A,B,V,D,es,es))
    T=MF.cubic_coefficients(A,B,V,D,F,E,Y)
    @test T*Nemo.matrix(F,10,1,m3)==Nemo.matrix(F,p*length(D),1,2qbilinear(A,B,V,D,es,ym)+cubic(A,D,es))
end

@testset "Core rank and jet results" begin
    r=kernel_report(5)
    @test (r["rank"],r["kernel"],r["gauge"],r["extra"])==(105,30,18,12)
    r,_,_,_=block_report(7)
    @test r["full_rank"]==271
    @test r["projected_block_dimensions"]==[25,41,41,41,41,41,41]
    F,A,B,V=fourier(5;kind=:eight_per_frequency,arithmetic=:numberfield)
    @test predict(A,B,V,allpairs(5))==target(5,allpairs(5),F)
    for n in (22,23)
        r=obstruction_report(n)
        @test r["training_rank"]==r["full_rank"]==105
        @test r["cubic_obstruction_rank_D"]==r["cubic_obstruction_rank_all"]==0
    end
end

@testset "Portable witness and exact existence certificate" begin
    path=joinpath(@__DIR__,"..","data","nongeneralizing_witness.txt")
    @test isfile(path)
    if isfile(path)
        x=load_witness(path);r=certify_interpolant(x)
        @test r["root_exists"] && r["rank_at_root"]==110
        @test r["heldout_error_lower_approx"]>0
        @test r["local_learning_coefficient"]==55
        mktemp() do name,io
            close(io);save_witness(name,x);@test load_witness(name)==x
        end
    end
end
