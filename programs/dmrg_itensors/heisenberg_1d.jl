#!/usr/bin/env julia
#
# 1D AFM Heisenberg model

include("common.jl")

function get_H(args)
    @assert args["L2"] == args["L"]
    @assert args["J22"] == 0
    @assert args["V"] == 0
    @assert args["h"] == 0

    L = args["L"]
    peri = args["peri"]
    mars = args["mars"]
    J2 = args["J2"]

    sites = siteinds("S=1/2", L; conserve_qns = args["zero_mag"])

    ampo = OpSum()

    function add_edge!(J, i1, i2)
        # Spin to Pauli
        J *= 4

        Jxy = 0.5 * J
        if mars
            d = i2 - i1
            if mod(d, 2) == 1
                Jxy *= -1
            end
        end

        i1 = mod1(i1, L)
        i2 = mod1(i2, L)
        ampo += J, "Sz", i1, "Sz", i2
        ampo += Jxy, "S+", i1, "S-", i2
        ampo += Jxy, "S-", i1, "S+", i2
    end

    for i = 1:L-1+peri
        add_edge!(1, i, i + 1)
    end

    if J2 != 0
        for i = 1:L-2+peri*2
            add_edge!(J2, i, i + 2)
        end
    end

    H = MPO(ampo, sites)

    return sites, H
end

function main()
    args = get_args("heis")
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
