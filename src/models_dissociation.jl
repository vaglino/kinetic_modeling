function cd3_dissociation!(du,u,p,t,dens,cons)
    # unpack rates and constants
    nᵣ,nₓ,n₃ = u
    k₋₃ = p[1]
    k₋₁,k₋₂ = cons
    # mᵣ,mₗ,mₓ,A = dens

    # model
    du[1] = dnᵣ = - k₋₁*nᵣ
    du[2] = dnₓ = - k₋₂*nₓ
    du[3] = dn₃ = - k₋₃*n₃

end

function bimolecular_diss_bell_activated_state_w_cons!(du,u,p,t,dens,cons;f=0.0)
    # unpack rates and constants
    Pᵣ,Pₐ = u
    x₋ᵒf,
        k₋ᵒs,x₋ᵒs,
        kᵒa,xᵒa,
        k₋ᵒa,x₋ᵒa = p
    k₋ᵒf = cons[1] # zero force diss, generally obtained from thermal fluct.

    k₋f = bell_diss_kT(k₋ᵒf, x₋ᵒf, f) # fast diss form L-R
    k₋s = bell_diss_kT(k₋ᵒs, x₋ᵒs, f) # slow diss from L-R*
    ka  = bell_diss_kT(kᵒa, xᵒa, f)   # activation from L-R to L-R*
    k₋a = bell_ass_kT(k₋ᵒa, x₋ᵒa, f)  # deactivation from L-R* to L-R

    # model
    du[1] = dPᵣ = - k₋f.*Pᵣ - ka.*Pᵣ + k₋a.*Pₐ
    du[2] = dPₐ = + ka.*Pᵣ - k₋a.*Pₐ - k₋s.*Pₐ
end

function cd3_tri_dissociation_bell!(du,u,p,t,dens,cons;f=0.0)
    # unpack rates and constants
    Pᵣ₁,Pₐ₁,
        Pᵣ₂,Pₐ₂,
        Pᵣ₃,Pₐ₃ = u

    x₋₃ᵒf,
        k₋₃ᵒs, x₋₃ᵒs,
        k₃ᵒa, x₃ᵒa,
        k₋₃ᵒa, x₋₃ᵒa = p

    #  for now declare the kᵢ from bimolecular as global variables
    k₋₁ᵒf, x₋₁ᵒf,
        k₋₁ᵒs, x₋₁ᵒs,
        k₁ᵒa, x₁ᵒa,
        k₋₁ᵒa, x₋₁ᵒa,
        k₋₂ᵒf, x₋₂ᵒf,
        k₋₂ᵒs, x₋₂ᵒs,
        k₂ᵒa, x₂ᵒa,
        k₋₂ᵒa, x₋₂ᵒa,
        k₋₃ᵒf         = cons



    k₋₁f = bell_diss_kT(k₋₁ᵒf, x₋₁ᵒf, f) # fast diss form L-R
    k₋₁s = bell_diss_kT(k₋₁ᵒs, x₋₁ᵒs, f) # slow diss from L-R*
    k₁a  = bell_diss_kT(k₁ᵒa, x₁ᵒa, f)   # activation from L-R to L-R*
    k₋₁a = bell_ass_kT(k₋₁ᵒa, x₋₁ᵒa, f)  # deactivation from L-R* to L-R

    k₋₂f = bell_diss_kT(k₋₂ᵒf, x₋₂ᵒf, f) # fast diss form L-R
    k₋₂s = bell_diss_kT(k₋₂ᵒs, x₋₂ᵒs, f) # slow diss from L-R*
    k₂a  = bell_diss_kT(k₂ᵒa, x₂ᵒa, f)   # activation from L-R to L-R*
    k₋₂a = bell_ass_kT(k₋₂ᵒa, x₋₂ᵒa, f)  # deactivation from L-R* to L-R

    k₋₃f = bell_diss_kT(k₋₃ᵒf, x₋₃ᵒf, f) # fast diss form L-R
    k₋₃s = bell_diss_kT(k₋₃ᵒs, x₋₃ᵒs, f) # slow diss from L-R*
    k₃a  = bell_diss_kT(k₃ᵒa, x₃ᵒa, f)   # activation from L-R to L-R*
    k₋₃a = bell_ass_kT(k₋₃ᵒa, x₋₃ᵒa, f)  # deactivation from L-R* to L-R

    # model
    # γϵ
    du[1] = dPᵣ₁ = - k₋₁f*Pᵣ₁ - k₁a*Pᵣ₁ + k₋₁a*Pₐ₁
    du[2] = dPₐ₁ = + k₁a*Pᵣ₁ - k₋₁a*Pₐ₁ - k₋₁s*Pₐ₁
    # δϵ
    du[3] = dPᵣ₂ = - k₋₂f*Pᵣ₂ - k₂a*Pᵣ₂ + k₋₂a*Pₐ₂
    du[4] = dPₐ₂ = + k₂a*Pᵣ₂ - k₋₂a*Pₐ₂ - k₋₂s*Pₐ₂
    # tri
    du[5] = dPᵣ₃ = - k₋₃f*Pᵣ₃ - k₃a*Pᵣ₃ + k₋₃a*Pₐ₃
    du[6] = dPₐ₃ = + k₃a*Pᵣ₃ - k₋₃a*Pₐ₃ - k₋₃s*Pₐ₃

end


## functions to calculate initial bond distribution vector (at t=0)

#=
if the delay between binding and beginning of the measurement is sufficiently
long, the probability of occupying either state will reach an equilibrium. This
equilibrium is independent of force
=#
function initial_state(kᵒa,k₋ᵒa)
        u₀1 = k₋ᵒa / (kᵒa + k₋ᵒa)
        u₀ = [u₀1, (1.0 - u₀1)]
end

#=
assumes that the time between binding and measurement is short relative to the
transition rates between the two bound states, then the bond will stay in the
initial state it occupied upon binding. The probabilities of occupying each of
these states can be derived from the principle of detailed balance.
=#
function initial_state(k₋ᵒf,k₋ᵒs,kᵒa,k₋ᵒa)
        u₀1 = (k₋ᵒa * k₋ᵒf) / (k₋ᵒa * k₋ᵒf + kᵒa * k₋ᵒs)
        u₀ = [u₀1, (1.0 - u₀1)]
end

#=
assumes enough time to reach equilibrium, and that the bond experiences some
fraction a of the measured force.
IMPORTANT: must add a*f to function call: initial_state(p, a * f)
=#


function initial_state(p::Vector,f)
        k₋ᵒf,x₋ᵒf,
        k₋ᵒs,x₋ᵒs,
        kᵒa,xᵒa,
        k₋ᵒa,x₋ᵒa = p

        k₋f = bell_diss_kT(k₋ᵒf, x₋ᵒf, f) # fast diss form L-R
        k₋s = bell_diss_kT(k₋ᵒs, x₋ᵒs, f) # slow diss from L-R*
        ka  = bell_diss_kT(kᵒa, xᵒa, f)   # activation from L-R to L-R*
        k₋a = bell_ass_kT(k₋ᵒa, x₋ᵒa, f)

        u₀1 = k₋a / (ka + k₋a)
        u₀ = [u₀1, (1.0 - u₀1)]
end

# function
# u = []
# for f in Fs
#         u₀ = initial_state(rates,f)
#         push!(u,u₀)
# end
#
# plot(reduce(hcat,u))

function initial_u(k,M)
    if ~isempty(M.cons)
        k = vcat(M.cons,k)
    end

    if M.u₀_type == "equilibrium"
        k₋ᵒf,x₋ᵒf,k₋ᵒs,x₋ᵒs,kᵒa,xᵒa,k₋ᵒa,x₋ᵒa = k
        new_u₀ = initial_state(kᵒa,k₋ᵒa)

    elseif M.u₀_type == "detailed_balance"
        k₋ᵒf,x₋ᵒf,k₋ᵒs,x₋ᵒs,kᵒa,xᵒa,k₋ᵒa,x₋ᵒa = k
        new_u₀ = initial_state(k₋ᵒf,k₋ᵒs,kᵒa,k₋ᵒa)

    elseif M.u₀_type == "force"
        f = M.cons
        a = k[end]
        new_u₀ = initial_state(k,a*f)

    else
        new_u₀ = []
    end

    if isempty(M.u₀)
        u₀ = new_u₀
    else
        u₀ = M.u₀
        u₀ = convert.(eltype(new_u₀),u₀)
        u₀ = append!(u₀, new_u₀) # append new initial states to provided initial states
    end
    return u₀
end
