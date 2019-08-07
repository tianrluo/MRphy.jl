using Test
using Unitful, UnitfulMR

using MRphy

@time begin

@testset "utils tests" for _ = [1]
  dim, k = (2,2), [1,2,3,4,0]u"cm^-1"
  gT, gR = k2g(k, true), k2g(k, false)
  @test CartesianLocations(dim) ==
        CartesianLocations(dim,false) .- collect(Tuple(ctrSub(dim)))' ==
        [-1 -1; 0 -1; -1 0; 0 0];
  @test ctrSub((1,2,3,4,5,6)) == CartesianIndex(1,2,2,3,3,4)
  @test ctrInd.([(3,4), (4,3)]) == [8, 7] # center index of fftshift
  @test g2k(gT, true) ≈ k && k2g(k, false) ≈ gR # `≈` in case of numeric errors
  @test g2s(gT, dt=1u"s") ≈ [gT[1]; diff(gT, dims=1)]/1u"s"
end

#= Core Features Tests =#
M, nt = [1 0 0; 0 1 0; 0 0 1], 512
t = 0:nt-1

rf = (10*(cos.(t/nt*2π) + sin.(t/nt*2π)im))u"Gauss"
gr = [zeros(nt,2) (10*atan.(t.-round(nt/2))/π)]u"Gauss/cm"
dt, des = 4e-6u"s", "test pulse"

p = Pulse(rf, gr; dt=dt, des=des)

spa = mSpinArray(trues((3,1)); M=M)

fov = Quantity.((3,3,1), u"cm")
cube = mSpinCube(trues((3,3,1)), fov)

@testset "mObjects tests" for _ = [1] # setting _ = [] shuts down this testset

  @testset "`AbstractPulse` constructor" for _ = [1]
    @test isequal(p, Pulse(p.rf, p.gr; dt=dt,des=des))
    @test p != Pulse(p.rf, p.gr; dt=dt,des=des)
  end

  @testset "`AbstractSpinArray` constructor" for _ = [1]
    @test isequal(spa, mSpinArray(spa.mask; M=spa.M))
    @test isequal(spa, mSpinArray(size(spa); M=spa.M))
  end

  @testset "`AbstractSpnCube` constructor" for _ = [1]
    @test isequal(cube, mSpinCube(cube.mask, fov))
    @test isequal(cube, mSpinCube(cube.dim, fov))
  end

end

spa.T1, spa.T2 = 1u"s", 4e-2u"s"
loc = [zeros((3,2)) ones((3,1))]u"cm"

Mo, _ = blochSim(spa, p, loc; doHist=false)

@testset "blochSim tests" for _ = [1] # setting _ = [] shuts down this testset

  B = cat(Pulse2B(p, loc)..., dims=3)
  @test B == cat(Pulse2B(p.rf, p.gr, loc)..., dims=3)

  _, Mhst = blochSim(M, B; T1=spa.T1,T2=spa.T2,γ=spa.γ,dt=p.dt,doHist=true);

  @test Mo            ≈
        Mhst[:,:,end] ≈
        [ 0.559535641648385  0.663342640621335  0.416341441715101;
          0.391994737048090  0.210182892388552 -0.860954821972489;
         -0.677062008711222  0.673391604920576 -0.143262993311057]

end

end

