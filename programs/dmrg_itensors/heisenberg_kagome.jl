#!/usr/bin/env julia
#
# Heisenberg model on kagome lattice with 18 spins

include("common.jl")

function get_H(args)
    @assert !args["peri"]
    @assert !args["mars"]
    @assert args["L"] == 0
    @assert args["L2"] == 0
    @assert args["J2"] == 0
    @assert args["J22"] == 0
    @assert args["V"] == 0
    @assert args["h"] == 0

    N = 18
    edges = [
        (0, 1),
        (0, 2),
        (0, 4),
        (0, 17),
        (1, 2),
        (1, 3),
        (1, 14),
        (2, 7),
        (2, 9),
        (3, 4),
        (3, 5),
        (3, 14),
        (4, 5),
        (4, 17),
        (5, 6),
        (5, 10),
        (6, 7),
        (6, 8),
        (6, 10),
        (7, 8),
        (7, 9),
        (8, 13),
        (8, 15),
        (9, 10),
        (9, 11),
        (10, 11),
        (11, 12),
        (11, 16),
        (12, 13),
        (12, 14),
        (12, 16),
        (13, 14),
        (13, 15),
        (15, 16),
        (15, 17),
        (16, 17),
    ]

    sites = siteinds("S=1/2", N; conserve_qns = args["zero_mag"])

    ampo = OpSum()

    function add_edge!(J, i, j)
        # Spin to Pauli
        J *= 4

        Jxy = 0.5 * J
        ampo += J, "Sz", i, "Sz", j
        ampo += Jxy, "S+", i, "S-", j
        ampo += Jxy, "S-", i, "S+", j
    end

    for (i, j) in edges
        add_edge!(1, i + 1, j + 1)
    end

    H = MPO(ampo, sites)

    return sites, H
end

function main()
    args = get_args("heis_kag")
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
