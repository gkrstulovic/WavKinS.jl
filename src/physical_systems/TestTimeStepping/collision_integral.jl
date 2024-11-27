using Base.Threads
# Here we define collisional integrals for TestTimeStepping system

@doc raw"""
    St_k!(Run::TestTimeStepping)

Compute collision integral for [`TestTimeStepping`](@ref "TestTimeStepping").

"""
function St_k!(Run::TestTimeStepping)
    Nk = Run.Nk
    nk = Nk.nk
    Sk = Run.Sk.nk

    @sync for id = 1:Threads.nthreads()
        inds = Run.partition[id]
        @spawn begin
            for i in inds
                Sk[i] = nk[i]^2
            end
        end
    end
end
