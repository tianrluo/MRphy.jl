using Test
using Unitful, UnitfulMR

using MRphy
using MRphy: NoUnitChk

@time begin

M, nt = [1 0 0; 0 1 0; 0 0 1], 512
t = 0:nt-1

rf = (10*(cos.(t/nt*2π) + sin.(t/nt*2π)im))u"Gauss"
gr = [zeros(nt,2) (10*atan.(t.-round(nt/2))/π)]u"Gauss/cm"
dt, des = 4e-6u"s", "test pulse"

p = Pulse(rf, gr; dt=dt, des=des)

spa = mSpinArray(trues((3,1)); M=M)

fov = Quantity.((3,3,1), u"cm")
cube = mSpinCube(trues((3,3,1)), fov)

@testset "mObjects tests" for i = [1] # setting i = [] shuts down this testset

  rf_ul, gr_ul, dt_ul = # _ul: unitless
  ustrip.(u"Gauss", rf), ustrip.(u"Gauss/cm", gr), ustrip.(u"s", dt)

  @testset "`AbstractPulse` constructor" for i = [1]
    @test isequal(p, Pulse(NoUnitChk(),rf_ul,gr_ul; dt=dt_ul,des=des))
    @test p != Pulse(NoUnitChk(),rf_ul,gr_ul; dt=dt_ul,des=des)
  end

  @testset "`AbstractSpinArray` constructor" for i = [1]
    @test isequal(spa, mSpinArray(NoUnitChk(), spa.mask; M=spa.M))
    @test isequal(spa, mSpinArray(size(spa); M=spa.M))
  end

  @testset "`AbstractSpnCube` constructor" for i = [1]
    @test isequal(cube, mSpinCube(NoUnitChk(), cube.mask, ustrip.(fov)))
    @test isequal(cube, mSpinCube(cube.dim, fov))
  end

end

spa.T1, spa.T2 = 1u"s", 4e-2u"s"
loc = [zeros((3,2)) ones((3,1))]u"cm"

Mo, _ = blochSim(spa, p, loc; doHist=false)

@testset "blochSim tests" for i = [1] # setting i = [] shuts down this testset

  B = Pulse2B(p, loc)
  @test B == Pulse2B(rf, gr, loc)

  _, Mhst = blochSim(M, B; T1=spa.T1,T2=spa.T2,γ=spa.γ,dt=p.dt,doHist=true);

  @test Mo            ≈
        Mhst[:,:,end] ≈
        [ 0.559535641648385  0.663342640621335  0.416341441715101;
          0.391994737048090  0.210182892388552 -0.860954821972489;
         -0.677062008711222  0.673391604920576 -0.143262993311057]

end

end

