
function trimolecular_1path_diff!(du,u,p,t,dens,cons)
    # unpack rates and constants
    nᵣ,nᵣₓ,n₃ = u
    kₓ,k₋ₓ,k₃,k₋₃ = p
    k₁,k₋₁ = cons
    mᵣ,mₗ,mₓ,A = dens

    # model
    du[1] = dnᵣ  = A*k₁*mᵣ*mₗ - k₋₁*nᵣ - A*kₓ*mₓ*nᵣ + k₋ₓ*nᵣₓ
    du[2] = dnᵣₓ = A*kₓ*mₓ*nᵣ - k₋ₓ*nᵣₓ - k₃*nᵣₓ + k₋₃*n₃
    du[3] = dn₃  = k₃*nᵣₓ - k₋₃*n₃

end
