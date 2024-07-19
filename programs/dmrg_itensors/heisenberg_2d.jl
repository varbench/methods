#!/usr/bin/env julia
#
# 2D AFM Heisenberg model

include("common.jl")

function get_H(args)
    @assert args["V"] == 0
    @assert args["h"] == 0

    L = args["L"]
    L2 = args["L2"]
    peri = args["peri"]
    mars = args["mars"]
    J2 = args["J2"]
    J22 = args["J22"]

    sites = siteinds("S=1/2", L * L2; conserve_qns = args["zero_mag"])

    ampo = OpSum()

    function ind(i, j)
        i = mod1(i, L)
        j = mod1(j, L2)
        # Snake ordering
        if i % 2 == 1
            return (i - 1) * L2 + j
        else
            return (i - 1) * L2 + (L2 + 1 - j)
        end
    end

    function add_edge!(J, i1, j1, i2, j2)
        # Spin to Pauli
        J *= 4

        Jxy = 0.5 * J
        if mars
            d = i2 - i1 + j2 - j1
            if mod(d, 2) == 1
                Jxy *= -1
            end
        end

        k1 = ind(i1, j1)
        k2 = ind(i2, j2)
        ampo += J, "Sz", k1, "Sz", k2
        ampo += Jxy, "S+", k1, "S-", k2
        ampo += Jxy, "S-", k1, "S+", k2
    end

    for i = 1:L
        for j = 1:L2-1+peri
            add_edge!(1, i, j, i, j + 1)
        end
    end
    for i = 1:L-1+peri
        for j = 1:L2
            add_edge!(1, i, j, i + 1, j)
        end
    end

    # Diagonal term, needed for triangular lattice and J1-J2 model
    if J2 != 0
        for i = 1:L-1+peri
            for j = 1:L2-1+peri
                add_edge!(J2, i, j, i + 1, j + 1)
            end
        end
    end

    # Inverse diagonal term, needed for J1-J2 model
    if J22 != 0
        for i = 1:L-1+peri
            for j = 1:L2-1+peri
                add_edge!(J22, i, j + 1, i + 1, j)
            end
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
