#!/usr/bin/env julia
#
# 1D AFM Ising model
# No zero magnetization constraint

include("common.jl")

function get_H(args)
    @assert args["L2"] == args["L"]
    @assert !args["mars"]
    @assert args["J2"] == 0
    @assert args["J22"] == 0
    @assert args["V"] == 0

    L = args["L"]
    peri = args["peri"]
    h = args["h"]

    @assert !args["zero_mag"]
    sites = siteinds("S=1/2", L)

    ampo = OpSum()

    function add_edge!(J, i1, i2)
        i1 = mod1(i1, L)
        i2 = mod1(i2, L)
        # Spin to Pauli
        ampo += 4 * J, "Sz", i1, "Sz", i2
    end

    for i = 1:L-1+peri
        add_edge!(1, i, i + 1)
    end

    if h != 0
        for i = 1:L
            # Spin to Pauli
            ampo += 2 * h, "Sx", i
        end
    end

    H = MPO(ampo, sites)

    return sites, H
end

function main()
    args = get_args("ising")
    @show args["out_filename"]

    sites, H = get_H(args)
    psi = get_init_psi_spin(args, sites)
    sweeps = get_sweeps(args)
    energy, psi = dmrg(H, psi, sweeps)
    @show energy

    h5open(args["out_filename"], "w") do f
        write(f, "psi", dense(psi))
    end

    GC.gc()

    H_sqr = real(inner(H, psi, H, psi))
    energy_var = H_sqr - energy^2
    @show energy_var
end

main()
