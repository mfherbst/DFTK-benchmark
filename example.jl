include("systems.jl")
include("payload.jl")

# Some setup wrapper of DFTK ...
setup_threading()

# In typical computational settings larger values for Ecut would be used.

# Iron and aluminium are "classical" examples.
# system = make_iron()
system = make_aluminium()

# This is a less representative for the plane-wave DFT use case
# system = make_caffeine()

# Print system info ...
if mpi_master()
    display("text/plain", system.basis)
end
guess = reproducible_guess(system)

# Compile code
payload(system.basis; guess..., maxiter=1)

# Run and time it
DFTK.reset_timer!(DFTK.timer)
payload(system.basis; guess..., maxiter=10)
println(DFTK.timer)
