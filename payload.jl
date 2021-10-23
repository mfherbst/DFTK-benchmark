"""
This sets up a consistent starting point for tests.
"""
function reproducible_guess(system::NamedTuple; tol=100)
    ρ0 = guess_density(system.basis, system.magnetic_moments)
    scfres = self_consistent_field(system.basis; ρ=ρ0, tol, callback=identity)
    (; ψ=scfres.ψ, ρ=scfres.ρ)
end


"""
The call to time and optimise.
"""
function payload(basis; ψ, ρ, tol=1e-14, kwargs...)
    mixargs = (; mixing=SimpleMixing()) # <- remove this Parameter to enable
    #                                        default settings in DFTK
    self_consistent_field(basis; ρ, ψ, tol, mixargs..., kwargs...)
end
