var documenterSearchIndex = {"docs":
[{"location":"#MRphy.jl-Documentation-1","page":"MRphy.jl Documentation","title":"MRphy.jl Documentation","text":"","category":"section"},{"location":"#Modules-1","page":"MRphy.jl Documentation","title":"Modules","text":"","category":"section"},{"location":"#","page":"MRphy.jl Documentation","title":"MRphy.jl Documentation","text":"MRphy","category":"page"},{"location":"#MRphy","page":"MRphy.jl Documentation","title":"MRphy","text":"General Comments:\n\nnM, number of spins, as magnetic spin vectors are often denoted as 𝑀.\nnT, number of steps/time-points.\n\nUnitful.jl related:\n\n𝐁 = 𝐌*𝐈^-1*𝐓^-2, dimension of magnetic field strength.\n𝐅 = 𝐓^-1, dimension of temporal frequency.\n𝐊 = 𝐋^-1, dimension of spatial frequency.\n𝚪 = 𝐅/𝐁, dimension of gyro ratio.\n\n\n\n\n\n","category":"module"},{"location":"#","page":"MRphy.jl Documentation","title":"MRphy.jl Documentation","text":"Modules = [MRphy.utils]\nOrder = [:module]","category":"page"},{"location":"#MRphy.utils","page":"MRphy.jl Documentation","title":"MRphy.utils","text":"Some utilities functions routinely used in MR simulations.\n\n\n\n\n\n","category":"module"},{"location":"#","page":"MRphy.jl Documentation","title":"MRphy.jl Documentation","text":"See Utilities","category":"page"},{"location":"#","page":"MRphy.jl Documentation","title":"MRphy.jl Documentation","text":"Modules = [MRphy.SteadyStates,\n           MRphy.SteadyStates.Signal,\n           MRphy.SteadyStates.RFSpoiling]\nOrder = [:module]","category":"page"},{"location":"#MRphy.SteadyStates","page":"MRphy.jl Documentation","title":"MRphy.SteadyStates","text":"Some steady state properties of common sequences.\n\n\n\n\n\n","category":"module"},{"location":"#MRphy.SteadyStates.Signal","page":"MRphy.jl Documentation","title":"MRphy.SteadyStates.Signal","text":"Analytical expressions of common steady states sequences signals.\n\n\n\n\n\n","category":"module"},{"location":"#MRphy.SteadyStates.RFSpoiling","page":"MRphy.jl Documentation","title":"MRphy.SteadyStates.RFSpoiling","text":"Tools for simulating RF spoiling in gradient echo sequences.\n\n\n\n\n\n","category":"module"},{"location":"#","page":"MRphy.jl Documentation","title":"MRphy.jl Documentation","text":"See SteadyStates","category":"page"},{"location":"#Constants-1","page":"MRphy.jl Documentation","title":"Constants","text":"","category":"section"},{"location":"#","page":"MRphy.jl Documentation","title":"MRphy.jl Documentation","text":"Modules = [MRphy]\nOrder = [:constant]","category":"page"},{"location":"#MRphy.γ¹H","page":"MRphy.jl Documentation","title":"MRphy.γ¹H","text":"const γ¹H = 4257.6u\"Hz/Gauss\"\n\nGyromagnetic ratio of water proton.\n\n\n\n\n\n","category":"constant"},{"location":"#Types-1","page":"MRphy.jl Documentation","title":"Types","text":"","category":"section"},{"location":"#","page":"MRphy.jl Documentation","title":"MRphy.jl Documentation","text":"Modules = [MRphy]\nOrder = [:type]","category":"page"},{"location":"#MRphy.AbstractPulse","page":"MRphy.jl Documentation","title":"MRphy.AbstractPulse","text":"An abstract type for pulses.\n\nSee also: Pulse.\n\n\n\n\n\n","category":"type"},{"location":"#MRphy.AbstractSpinArray","page":"MRphy.jl Documentation","title":"MRphy.AbstractSpinArray","text":"This type keeps the essentials of magnetic spins. Its instance struct must contain all fields listed listed in the exemplary struct mSpinArray.\n\nMisc\n\nMight make AbstractSpinArray <: AbstractArray in a future version\n\nSee also: mSpinArray, AbstractSpinCube.\n\n\n\n\n\n","category":"type"},{"location":"#MRphy.AbstractSpinCube","page":"MRphy.jl Documentation","title":"MRphy.AbstractSpinCube","text":"AbstractSpinCube <: AbstractSpinArray. This type inherits AbstractSpinArray as a field. Its instance struct must contain all fields listed in the exemplary struct mSpinCube.\n\nSee also: AbstractSpinArray, mSpinCube.\n\n\n\n\n\n","category":"type"},{"location":"#MRphy.B0D","page":"MRphy.jl Documentation","title":"MRphy.B0D","text":"B0D = Quantity{<:Real, 𝐁}\n\nType of magetic field strength. Based on Unitful.Quantity.\n\nExamples:\n\njulia> (1u\"Gauss\")::B0D\n1 Gauss\n\n\n\n\n\n","category":"type"},{"location":"#MRphy.F0D","page":"MRphy.jl Documentation","title":"MRphy.F0D","text":"F0D =  Quantity{<:Real, 𝐅}\n\nType of temporal frequency. Based on Unitful.Quantity.\n\nExamples:\n\njulia> (1u\"s^-1\")::F0D\n1 s^-1\n\n\n\n\n\n","category":"type"},{"location":"#MRphy.GR0D","page":"MRphy.jl Documentation","title":"MRphy.GR0D","text":"GR0D = Quantity{<:Real, 𝐁/𝐋}\n\nType of magnetic gradient. Based on Unitful.Quantity.\n\nExamples:\n\njulia> (1u\"Gauss/cm\")::GR0D\n1 Gauss cm^-1\n\n\n\n\n\n","category":"type"},{"location":"#MRphy.K0D","page":"MRphy.jl Documentation","title":"MRphy.K0D","text":"K0D =  Quantity{<:Real, 𝐊}\n\nType of spatial frequency. Based on Unitful.Quantity.\n\nExamples:\n\njulia> (1u\"cm^-1\")::K0D\n1 cm^-1\n\n\n\n\n\n","category":"type"},{"location":"#MRphy.L0D","page":"MRphy.jl Documentation","title":"MRphy.L0D","text":"L0D = Quantity{<:Real, 𝐋}\n\nType of length. Based on Unitful.Quantity.\n\nExamples:\n\njulia> (1u\"cm\")::L0D\n1 cm\n\n\n\n\n\n","category":"type"},{"location":"#MRphy.Pulse","page":"MRphy.jl Documentation","title":"MRphy.Pulse","text":"Pulse(rf, gr; dt=(4e-6)u\"s\", des=\"generic pulse\")\n\nA struct for typical MR pulses: Pulse <: AbstractPulse.\n\nFields:\n\nMutable:\n\nrf::TypeND(RF0D, [1,2]) (nT,) or (nT, nCoils).\ngr::TypeND(GR0D, [2]) (nT, 3), where 3 accounts for x-y-z channels.\ndt::T0D (1,), simulation temporal step size, i.e., dwell time.\ndes::String, an description of the pulse to be constructed.\n\nSee also: AbstractPulse.\n\n\n\n\n\n","category":"type"},{"location":"#MRphy.RF0D","page":"MRphy.jl Documentation","title":"MRphy.RF0D","text":"RF0D = Quantity{<:Union{Real, Complex}, 𝐁}\n\nType of magnetic RF. Based on Unitful.Quantity.\n\nExamples:\n\njulia> ((1+1im)u\"Gauss\")::RF0D\n1 + 1im Gauss\n\n\n\n\n\n","category":"type"},{"location":"#MRphy.T0D","page":"MRphy.jl Documentation","title":"MRphy.T0D","text":"T0D = Quantity{<:Real, 𝐓}\n\nType of time. Based on Unitful.Quantity.\n\nExamples:\n\njulia> (1u\"s\")::T0D\n1 s\n\n\n\n\n\n","category":"type"},{"location":"#MRphy.mSpinArray","page":"MRphy.jl Documentation","title":"MRphy.mSpinArray","text":"mSpinArray(dim::Dims; T1=1.47u\"s\", T2=0.07u\"s\", γ=γ¹H, M=[0. 0. 1.])\nmSpinArray(mask::BitArray; T1=1.47u\"s\", T2=0.07u\"s\", γ=γ¹H, M=[0. 0. 1.])\n\nAn exemplary struct instantiating AbstractSpinArray.\n\nFields:\n\nImmutable:\n\ndim::Dims (nd,): nM ← prod(dim), dimension of the object.\nmask::BitArray (nx,(ny,(nz))): Mask for M, dim == (nx,(ny,(nz)))\n\nMutable:\n\nT1::TypeND(T0D, [0,1]) (1,) or (nM,): Longitudinal relaxation coeff.\nT2::TypeND(T0D, [0,1]) (1,) or (nM,): Transversal relaxation coeff.\nγ::TypeND(Γ0D, [0,1])  (1,) or (nM,): Gyromagnetic ratio.\nM::TypeND(Real, [2])   (count(mask), 3):  Magnetic spins, (𝑀x,𝑀y,𝑀z).\n\nNotes:\n\noff-resonance, Δf, and locations, loc, are intentionally unincluded, as they are not intrinsic to spins, and can change over time. Unincluding them allows extensional subtypes specialized for, e.g., arterial spin labelling.\n\nSee also: AbstractSpinArray.\n\n\n\n\n\n","category":"type"},{"location":"#MRphy.mSpinCube","page":"MRphy.jl Documentation","title":"MRphy.mSpinCube","text":"spincube = mSpinCube(dim::Dims{3}, fov; ofst, Δf, T1, T2, γ)\nspincube = mSpinCube(mask::BitArray{3}, fov; ofst, Δf, T1, T2, γ)\n\ndim, mask, T1, T2, and γ are passed to mSpinArray constructors.\n\nAn exemplary struct instantiating AbstractSpinCube, designed to model a set of regularly spaced spins, e.g., a volume.\n\nFields:\n\nImmutable:\n\nspinarray::AbstractSpinArray (1,): inherited AbstractSpinArray struct\nfov::NTuple{3,L0D} (3,): field of view.\nofst::NTuple{3,L0D} (3,): fov offset from magnetic field iso-center.\nloc::TypeND(L0D, [2])  (nM, 3): location of spins.\n\nMutable:\n\nΔf::TypeND(F0D, [0,1]) (1,) or (nM,): off-resonance map.\n\nSee also: AbstractSpinCube.\n\n\n\n\n\n","category":"type"},{"location":"#MRphy.Γ0D","page":"MRphy.jl Documentation","title":"MRphy.Γ0D","text":"Γ0D = Quantity{<:Real, 𝚪}\n\nType of gyro magnetic ratio. Based on Unitful.Quantity.\n\nExamples:\n\njulia> (1u\"Hz/Gauss\")::Γ0D\n1 Hz Gauss^-1\n\n\n\n\n\n","category":"type"},{"location":"#MRphy.AbstractSpinBolus","page":"MRphy.jl Documentation","title":"MRphy.AbstractSpinBolus","text":"UNDER CONSTRUCTION\n\nAbstractSpinBolus <: AbstractSpinArray. This type inherits AbstractSpinArray as a field. Its instance struct must contain all fields listed in the exemplary struct mSpinBolus.\n\nSee also: AbstractSpinArray, mSpinBolus.\n\n\n\n\n\n","category":"type"},{"location":"#MRphy.mSpinBolus","page":"MRphy.jl Documentation","title":"MRphy.mSpinBolus","text":"UNDER CONSTRUCTION\n\nAn exemplary struct instantiating AbstractSpinBolus, designed to model a set of moving spins, e.g., a blood bolus in ASL context.\n\nSee also: AbstractSpinBolus.\n\n\n\n\n\n","category":"type"},{"location":"#Functions-1","page":"MRphy.jl Documentation","title":"Functions","text":"","category":"section"},{"location":"#","page":"MRphy.jl Documentation","title":"MRphy.jl Documentation","text":"TypeND\nrfgr2B\nPulse2B\nblochSim\nblochSim!\nB2AB\nB2UΦ\nB2UΦ!\nUΦRot\nUΦRot!\napplyPulse\napplyPulse!\nfreePrec\nfreePrec!\nMRphy.ExceptionImmutableField","category":"page"},{"location":"#MRphy.TypeND","page":"MRphy.jl Documentation","title":"MRphy.TypeND","text":"TypeND(T,Ns) = Union{AbstractArray{<:T,Ns[1]}, AbstractArray{<:T,Ns[2]},...}\n\nSugar for creating Union{<:T typed array of different dimensions}.\n\nUsage\n\nINPUTS:\n\nT::Type (1,), the underlying type of the union.\nNs::Array{Int64,1} (# diff dims,), an array of wanted dimensions.\n\n\n\n\n\nTypeND(T, ::Colon) = AbstractArray{<:T}\n\nSugar for creating <:T typed array of arbitrary dimensions.\n\n\n\n\n\n","category":"function"},{"location":"#MRphy.rfgr2B","page":"MRphy.jl Documentation","title":"MRphy.rfgr2B","text":"B = rfgr2B(rf, gr, loc=[0 0 0]u\"cm\"; Δf=0u\"Hz\", b1Map=1, γ=γ¹H)\n\nTurn rf, rf, and gradient, gr, into 𝐵-effective magnetic field.\n\nINPUTS:\n\nrf::TypeND(RF0D, [1,2]) (nT, (nCoil))\ngr::TypeND(GR0D, [2])   (nT, 3)\nloc::TypeND(L0D, [2])   (1,3) or (nM, 3), locations.\n\nKEYWORDS:\n\nΔf::TypeND(F0D, [0,1,2]) (1,)  or (nM,), off-resonance.\nb1Map::TypeND(Union{Real,Complex},[0,1,2]) (1,) or (nM,(nCoils)),  transmit sensitivity.\nγ::TypeND(Γ0D, [0,1]) (1,)  or (nM,), gyro-ratio\n\nOUTPUS:\n\nB, generator of TypeND(B0D, [2]) (1,1,nT), 𝐵 field.\n\nSee also: Pulse2B, blochSim.\n\nTODO:\n\nSupport loc, Δf, and b1Map being Base.Generators.\n\n\n\n\n\n","category":"function"},{"location":"#MRphy.Pulse2B","page":"MRphy.jl Documentation","title":"MRphy.Pulse2B","text":"B = Pulse2B(pulse::Pulse, loc; Δf, b1Map, γ)\n\nCreate effective magnetic field, 𝐵, from input pulse.\n\nSee also: rfgr2B, B2UΦ, blochSim.\n\n\n\n\n\nB = Pulse2B(pulse::Pulse, spa::AbstractSpinArray, loc; Δf, b1Map)\n\n...with γ=spa.γ.\n\n\n\n\n\nB = Pulse2B(pulse::Pulse, cb::AbstractSpinCube; b1Map)\n\n...with loc, Δf, γ = cb.loc, cb.Δf, cb.γ.\n\n\n\n\n\n","category":"function"},{"location":"#MRphy.blochSim","page":"MRphy.jl Documentation","title":"MRphy.blochSim","text":"blochSim(M, B; T1, T2, γ, dt, doHist)\n\nSame as blochSim!(M, B; T1,T2,γ,dt,doHist), M will not be updated.\n\nSee also: blochSim\n\n\n\n\n\nblochSim(M, A, B)\n\nSame as blochSim(M, A, B), M will not be updated.\n\n\n\n\n\n","category":"function"},{"location":"#MRphy.blochSim!","page":"MRphy.jl Documentation","title":"MRphy.blochSim!","text":"blochSim!(M, B; T1=(Inf)u\"s\",T2=(Inf)u\"s\",γ=γ¹H,dt=(4e-6)u\"s\",doHist=false)\n\nOld school 𝐵-effective magnetic field, B, based bloch simulation. Globally or spin-wisely apply B over spins, M. M will be updated by the results.\n\nINPUTS:\n\nM::TypeND(Real, [2]) (nM, xyz): input spins' magnetizations.\nB::Union{TypeND(B0D, [2,3]), Base.Generator}: Global, (nT,xyz); Spin-wise, (nM,xyz,nT).\n\nKEYWORDS:\n\nT1 & T2 ::TypeND(T0D, [0,1]): Global, (1,); Spin-wise, (nM,1).\nγ::TypeND(Γ0D, [0,1]): Global, (1,); Spin-wise, (nM, 1). gyro ratio\ndt::T0D (1,), simulation temporal step size, i.e., dwell time.\ndoHist::Bool, whether to output spin history through out B.\n\nOUTPUTS:\n\nM::TypeND(Real, [2]) (nM, xyz): spins after applying B.\nMhst::TypeND(Real, [3]) (nM, xyz, nT): spins history during B.\n\nSee also: applyPulse, blochSim!.\n\nNotes:\n\nNot much sanity check inside this function, user is responsible for matching up the dimensions.\nPut relax at the end of each time step may still be inaccurate, since physically spins relax continuously, this noise/nuance may worth study for applications like fingerprinting simulations, etc.\n\n\n\n\n\nblochSim!(M, A, B)\n\nHargreave's 𝐴/𝐵, mat/vec, based bloch simulation. Globally or spin-wisely apply matrix A and vector B over spins, M, described in doi:10.1002/mrm.1170\n\nINPUTS:\n\nM::TypeND(Real, [2]) (nM, xyz): input spins' magnetizations.\nA::TypeND(AbstractFloat,[3]) (nM, 3,3), A[iM,:,:] is the iM-th 𝐴.\nB::TypeND(AbstractFloat,[2]) (nM, 3), B[iM,:] is the iM-th 𝐵.\n\nOUTPUTS:\n\nM::TypeND(Real, [2]) (nM, xyz): spins after applying B.\n\n\n\n\n\n","category":"function"},{"location":"#MRphy.B2AB","page":"MRphy.jl Documentation","title":"MRphy.B2AB","text":"B2AB(B; T1=(Inf)u\"s\", T2=(Inf)u\"s\", γ=γ¹H, dt=(4e-6)u\"s\")\n\nTurn B-effective into Hargreave's 𝐴/𝐵, mat/vec, see: doi:10.1002/mrm.1170.\n\nINPUTS:\n\nB::Union{TypeND(B0D, [2,3]), Base.Generator}: Global, (nT,xyz); Spin-wise, (nM,xyz,nT).\n\nKEYWORDS:\n\nT1 & T2 ::TypeND(T0D, [0,1]): Global, (1,); Spin-wise, (nM,1).\nγ::TypeND(Γ0D, [0,1]): Global, (1,); Spin-wise, (nM, 1). gyro ratio\ndt::T0D (1,), simulation temporal step size, i.e., dwell time.\n\nOUTPUTS:\n\nA::TypeND(AbstractFloat,[3]) (nM, 3,3), A[iM,:,:] is the iM-th 𝐴.\nB::TypeND(AbstractFloat,[2]) (nM, 3), B[iM,:] is the iM-th 𝐵.\n\nSee also: rfgr2B, Pulse2B.\n\n\n\n\n\n","category":"function"},{"location":"#MRphy.B2UΦ","page":"MRphy.jl Documentation","title":"MRphy.B2UΦ","text":"B2UΦ(B::TypeND(B0D,[2,3]); γ::TypeND(Γ0D,[0,1]), dt::T0D=4e-6u\"s\")\n\nGiven 𝐵-effective, B, compute rotation axis/angle, U/Φ.\n\nINPUTS:\n\nB::TypeND(B0D, [2,3]) (1,3,nT) or (nM, 3, nT), 𝐵 field.\n\nKEYWORDS:\n\nγ::TypeND(Γ0D, [0,1]): Global, (1,); Spin-wise, (nM, 1). gyro ratio\ndt::T0D (1,), simulation temporal step size, i.e., dwell time.\n\nOUTPUTS:\n\nU::TypeND(Real, [2,3]) (1,3,(nT)) or (nM,3,(nT)), axis.\nΦ::TypeND(Real, [2,3]) (1,1,(nT)) or (nM,1,(nT)), angle.\n\nSee also: B2UΦ!, UΦRot!.\n\nNotes:\n\nSomehow, in-place version, B2UΦ!(B,U,Φ; γ,dt), provokes more allocs in julia.\n\n\n\n\n\n","category":"function"},{"location":"#MRphy.B2UΦ!","page":"MRphy.jl Documentation","title":"MRphy.B2UΦ!","text":"B2UΦ!(B, U; Φ, γ, dt=(4e-6)u\"s\")\n\nIn-place version of B2UΦ. Somehow, B2UΦ!, provokes more allocs in julia.\n\nSee also: B2UΦ, blochSim.\n\n\n\n\n\n","category":"function"},{"location":"#MRphy.UΦRot","page":"MRphy.jl Documentation","title":"MRphy.UΦRot","text":"UΦRot(U, Φ, V)\n\nSame as UΦRot!(U, Φ, V, R), except not in-place.\n\nSee also: UΦRot!.\n\n\n\n\n\n","category":"function"},{"location":"#MRphy.UΦRot!","page":"MRphy.jl Documentation","title":"MRphy.UΦRot!","text":"UΦRot!(U, Φ, V, R)\n\nApply axis-angle, U-Φ based rotation on V. Rotation is broadcasted on V along its 3rd dimension. Results will overwrite into R.\n\nINPUTS:\n\nU::TypeND(AbstractFloat,[2]) (nM, 3), rotation axes in 3D, assumed unitary;\nΦ::TypeND(AbstractFloat,[1]) (nM,), rotation angles;\nV::TypeND(AbstractFloat,[2,3]) (nM, 3, (3)), vectors to be rotated;\nR::TypeND(AbstractFloat,[2,3]) (nM, 3, (3)), vectors rotated, i.e., results;\n\nOUTPUTS:\n\nR the input container R is also returned for convenience.\n\nSee also: UΦRot, B2UΦ!.\n\n\n\n\n\n","category":"function"},{"location":"#MRphy.applyPulse","page":"MRphy.jl Documentation","title":"MRphy.applyPulse","text":"applyPulse(spa::AbstractSpinArray, p::Pulse, loc; Δf, b1Map, doHist)\n\nTurn p into 𝐵-effective and apply it on spa.M, using its own M, T1, T2, γ.\n\nSee also: blochSim, freePrec.\n\n\n\n\n\napplyPulse(cb::AbstractSpinCube, p::Pulse; b1Map, doHist)\n\nTurn p into 𝐵-effective and apply it on cb.M, using its own M, T1, T2, γ.\n\n\n\n\n\n","category":"function"},{"location":"#MRphy.applyPulse!","page":"MRphy.jl Documentation","title":"MRphy.applyPulse!","text":"applyPulse!(spa::AbstractSpinArray, p::Pulse, loc; Δf, b1Map, doHist)\n\nUpdate spa.M before return.\n\n\n\n\n\napplyPulse!(cb::AbstractSpinCube, p::Pulse; b1Map, doHist)\n\nUpdate cb.M before return.\n\n\n\n\n\n","category":"function"},{"location":"#MRphy.freePrec","page":"MRphy.jl Documentation","title":"MRphy.freePrec","text":"freePrec(M, t; Δf, T1, T2)\n\nSame as freePrec!(M, t; Δf, T1, T2), M will not be updated.\n\nSee also: freePrec!.\n\n\n\n\n\nfreePrec(spa::AbstractSpinArray, t; Δf)\n\nspa::AbstractSpinArray free precess by t. spa.M will not be updated.\n\nSee also: applyPulse, freePrec!.\n\n\n\n\n\nfreePrec(cb::AbstractSpinCube, t)\n\ncb::AbstractSpinCube free precess by t. cb.M will not be updated.\n\nSee also: applyPulse, freePrec.\n\n\n\n\n\n","category":"function"},{"location":"#MRphy.freePrec!","page":"MRphy.jl Documentation","title":"MRphy.freePrec!","text":"freePrec!(M, t; Δf=0u\"Hz\", T1=(Inf)u\"s\", T2=(Inf)u\"s\")\n\nSpins, M, free precess by time t. M will be updated by the results.\n\nINPUTS:\n\nM::TypeND(Real, [2]) (nM, xyz): input spins' magnetizations.\nt::T0D (1,): duration of free precession.\n\nKEYWORDS:\n\nT1 & T2 ::TypeND(T0D, [0,1]): Global, (1,); Spin-wise, (nM,1).\n\nOUTPUTS:\n\nM::TypeND(Real, [2]) (nM, xyz): output spins' magnetizations.\n\nSee also: freePrec.\n\n\n\n\n\nfreePrec!(spa::AbstractSpinArray, t; Δf)\n\n...spa.M will updated by the results.\n\n\n\n\n\nfreePrec!(cb::AbstractSpinCube, t)\n\n...cb.M will be updated by the results.\n\n\n\n\n\n","category":"function"},{"location":"#MRphy.ExceptionImmutableField","page":"MRphy.jl Documentation","title":"MRphy.ExceptionImmutableField","text":"Throw ArgumentError when $x is an immutable field.\n\n\n\n\n\n","category":"function"},{"location":"#Utilities-1","page":"MRphy.jl Documentation","title":"Utilities","text":"","category":"section"},{"location":"#","page":"MRphy.jl Documentation","title":"MRphy.jl Documentation","text":"Modules = [MRphy.utils]\nOrder = [:function]","category":"page"},{"location":"#MRphy.utils.CartesianLocations","page":"MRphy.jl Documentation","title":"MRphy.utils.CartesianLocations","text":"CartesianLocations(dim::Dims, doShift::Bool=true)\n\nRetuns a (prod(dim),length(dim)) sized array of grid locations. doShift shifts the locations to be consistent with fftshift.\n\nExamples\n\njulia> loc = [CartesianLocations((2,2)), CartesianLocations((2,2),false)]\n2-element Array{Array{Int64,2},1}:\n [-1 -1; 0 -1; -1 0; 0 0]\n [1 1; 2 1; 1 2; 2 2]\n\n\n\n\n\n","category":"function"},{"location":"#MRphy.utils.ctrInd-Tuple{Tuple{Vararg{Int64,N}} where N}","page":"MRphy.jl Documentation","title":"MRphy.utils.ctrInd","text":"ctrInd(dim::Dims) = sum((dim.÷2) .* [1; cumprod([dim[1:end-1]...])])+1\n\nAs a separate fn, ensure consistent behariour of getting the linear index to the center of a Nd-array of size dim. This center should match fftshift's center.\n\nSee also: ctrSub.\n\n\n\n\n\n","category":"method"},{"location":"#MRphy.utils.ctrSub-Tuple{Tuple{Vararg{Int64,N}} where N}","page":"MRphy.jl Documentation","title":"MRphy.utils.ctrSub","text":"ctrSub(dim::Dims) = CartesianIndex(dim .÷ 2 .+ 1)\n\nAs a separate function, ensure consistent behaviour of getting ::CartesianIndex to the center of a Nd-Array of size dim. This center should match fftshift's center.\n\nSee also: ctrInd.\n\nNotes:\n\nThe function may be removed once julia FFT packages provides this functionality.\n\n\n\n\n\n","category":"method"},{"location":"#MRphy.utils.g2k","page":"MRphy.jl Documentation","title":"MRphy.utils.g2k","text":"g2k(g::TypeND(GR0D,:); isTx::Bool=false, dt::T0D=4e-6u\"s\", γ::Γ0D=γ¹H)\n\nCompute k-space from gradient.\n\nUsage\n\nINPUTS:\n\ng::TypeND(GR0D, :) (nSteps, Nd...), gradient\nisTx::Bool, if true, compute transmit k-space, k, ends at the origin.\n\nKEYWORDS:\n\ndt::T0D (1,), gradient temporal step size, i.e., dwell time.\nγ::Γ0D (1,), gyro-ratio.\n\nOUTPUTS:\n\nk::TypeND(K0D, :) (nSteps, Nd...), k-space, w/ unit u\"cm^-1\".\n\nSee also: k2g, g2s.\n\n\n\n\n\n","category":"function"},{"location":"#MRphy.utils.g2s-Tuple{AbstractArray{#s12,N} where N where #s12<:(Unitful.Quantity{#s12,𝐌*𝐈^-1*𝐋^-1*𝐓^-2,U} where U where #s12<:Real)}","page":"MRphy.jl Documentation","title":"MRphy.utils.g2s","text":"g2s(g::TypeND(GR0D,:); dt::T0D=4e-6u\"s\")\n\nSlew rate sl, of the gradient, g.\n\nUsage\n\nINPUTS:\n\ng::TypeND(GR0D, :) (nSteps, Nd...)\n\nKEYWORDS:\n\ndt::T0D (1,), gradient temporal step size, i.e., dwell time.\n\nOUTPUTS:\n\nsl::TypeND(Quantity{<:Real, 𝐁/𝐋/𝐓, :) (nSteps, Nd...), slew rate\n\nSee also: g2k, k2g.\n\nNote\n\nNo s2g is provided for the moment.\n\n\n\n\n\n","category":"method"},{"location":"#MRphy.utils.k2g","page":"MRphy.jl Documentation","title":"MRphy.utils.k2g","text":"k2g(k::TypeND(K0D,:), isTx::Bool=false; dt::T0D=4e-6u\"s\", γ::Γ0D=γ¹H)\n\nGradient, g, of the TxRx k-space, (trasmit/receive, excitation/imaging).\n\nUsage\n\nINPUTS:\n\nk::TypeND(K0D, :) (nSteps, Nd...), Tx or Rx k-space, w/ unit u\"cm^-1\".\nisTx::Bool, if true, compute transmit k-space, k, ends at the origin.\n\nKEYWORDS:\n\ndt::T0D (1,), gradient temporal step size, i.e., dwell time.\nγ::Γ0D (1,), gyro-ratio.\n\nOUTPUTS:\n\ng::TypeND(GR0D, :) (nSteps, Nd...), gradient\n\nNote\n\nThe function asserts if k ends at the origin for isTx==true.\n\nSee also: g2k, g2s.\n\n\n\n\n\n","category":"function"},{"location":"#SteadyStates-1","page":"MRphy.jl Documentation","title":"SteadyStates","text":"","category":"section"},{"location":"#","page":"MRphy.jl Documentation","title":"MRphy.jl Documentation","text":"Modules = [MRphy.SteadyStates,\n           MRphy.SteadyStates.Signal,\n           MRphy.SteadyStates.RFSpoiling]\nOrder = [:function]","category":"page"},{"location":"#MRphy.SteadyStates.Signal.SPGR-Tuple{Real}","page":"MRphy.jl Documentation","title":"MRphy.SteadyStates.Signal.SPGR","text":"SPGR(α; TR, T1)\n\nTE=0 Steady state SPGR signal. 10.1002/mrm.1910130109, eq.(1), ideal spoiling.\n\nINPUTS:\n\nα::Real (1,), tip angle in degree;\n\nKEYWORDS:\n\nTR::T0D (1,), repetition time;\nT1::T0D (1,), longitudinal relaxation coefficient;\n\nOUTPUTS:\n\nsig::Real (1,), steady-state signal.\n\nSee also: bSSFP, STFR.\n\n\n\n\n\n","category":"method"},{"location":"#MRphy.SteadyStates.Signal.STFR-Tuple{Real,Real}","page":"MRphy.jl Documentation","title":"MRphy.SteadyStates.Signal.STFR","text":"STFR(α, β; ϕ, Δf, T1, T2, Tg, Tf)\n\nTE=0 Steady state STFR signal. 10.1002/mrm.25146, eq.(2), ideal spoiling.\n\nINPUTS:\n\nα::Real (1,), tip-down angle in degree;\nβ::Real (1,), tip-up angle in degree;\n\nKEYWORDS:\n\nϕ::Real (1,), phase of the tip-up pulse in radians;\nΔf::F0D (1,), off-resonance in Hz;\nT1::T0D (1,), longitudinal relaxation coefficient;\nT2::T0D (1,), transverse relaxation coefficient;\nTg::T0D (1,), duration of gradient crusher;\nTf::T0D (1,), duration of free precession in each TR;\n\nOUTPUTS:\n\nsig::Real (1,), steady-state signal.\n\nSee also: bSSFP, SPGR.\n\n\n\n\n\n","category":"method"},{"location":"#MRphy.SteadyStates.Signal.bSSFP-Tuple{Real}","page":"MRphy.jl Documentation","title":"MRphy.SteadyStates.Signal.bSSFP","text":"bSSFP(α; TR, Δf, T1, T2)\n\nTE=0 Steady state bSSFP signal. 10.1002/jmri.24163, eq.(4), with ϕ=2π*Δf*TR.\n\nINPUTS:\n\nα::Real (1,), tip angle in degree;\n\nKEYWORDS:\n\nTR::T0D (1,), repetition time;\nΔf::F0D (1,), off-resonance in Hz;\nT1::T0D (1,), longitudinal relaxation coefficient;\nT2::T0D (1,), transverse relaxation coefficient;\n\nOUTPUTS:\n\nsig::Complex (1,), steady-state signal.\n\nSee also: SPGR, STFR.\n\n\n\n\n\n","category":"method"},{"location":"#MRphy.SteadyStates.RFSpoiling.FZstates-Tuple{AbstractArray{D,2} where D<:Real,Real}","page":"MRphy.jl Documentation","title":"MRphy.SteadyStates.RFSpoiling.FZstates","text":"FZstates(Φ, α; TR, T1, T2, FZ)\n\n𝐹, 𝑍 from: 10.1002/(SICI)1099-0534(1999)11:5<291::AID-CMR2>3.0.CO;2-J, eq.(7,8). eq.(7) refined to ÷√(2), instead of ÷2, as in 10.1002/mrm20736: eq.(2).\n\nAssume constant gradient spoiling of m⋅2π dephasing in each TR, m∈ℤ. In practice, if dephase is a constant but not exactly m⋅2π, the resulting states can be computed by convolving a sinc with the m⋅2π dephased results.\n\nINPUTS:\n\nΦ::TypeND(Real,[2]) (nC,nTR), nC: #C as C in QuadPhase. Typically, one simulates a range of Cs, picking a C yielding a signal intensity equals to that of ideal spgr spoiling.\nα::Real (1,), flip-angle.\n\nKEYWORDS:\n\nTR::T0D (1,), repetition time;\nT1::T0D (1,), longitudinal relaxation coefficient;\nT2::T0D (1,), transverse relaxation coefficient;\nFZ::NamedTuple, (Fs,Fcs,Zs), simulate from prescribed states if given:\nFs ::TypeND(Complex,[2]), transversal dephasing states, 𝐹ₙ;\nFcs::TypeND(Complex,[2]), conjugate transversal dephasing states, 𝐹₋ₙ*;\nZs ::TypeND(Complex,[2]), longitudinal states, 𝑍ₙ;\n\nOUTPUTS:\n\nFZ::NamedTuple, (Fs,Fcs,Zs), simulation results.\n\n\n\n\n\n","category":"method"},{"location":"#MRphy.SteadyStates.RFSpoiling.QuadPhase","page":"MRphy.jl Documentation","title":"MRphy.SteadyStates.RFSpoiling.QuadPhase","text":"QuadPhase(nTR::Integer, C::Real, B::Real, A::Real)\n\nQuadratically cycling phases in (°): Φ(n) = mod.(C⋅n² + B⋅n + A, 360), n=0:nTR-1\n\n\n\n\n\n","category":"function"},{"location":"#Index-1","page":"MRphy.jl Documentation","title":"Index","text":"","category":"section"},{"location":"#","page":"MRphy.jl Documentation","title":"MRphy.jl Documentation","text":"","category":"page"}]
}
