#!/usr/bin/env julia
#
# 2D t-V model

include("common.jl")

function get_H(args)
    @assert !args["mars"]
    @assert args["J2"] == 0
    @assert args["J22"] == 0
    @assert args["h"] == 0

    L = args["L"]
    L2 = args["L2"]
    peri = args["peri"]
    V = args["V"]

    sites = siteinds("Fermion", L * L2; conserve_qns = true)

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

    function add_edge!(i1, j1, i2, j2)
        k1 = ind(i1, j1)
        k2 = ind(i2, j2)
        ampo -= 1, "Cdag", k1, "C", k2
        ampo += 1, "C", k1, "Cdag", k2
        ampo += V, "N", k1, "N", k2
    end

    for i = 1:L
        for j = 1:L2-1+peri
            add_edge!(i, j, i, j + 1)
        end
    end
    for i = 1:L-1+peri
        for j = 1:L2
            add_edge!(i, j, i + 1, j)
        end
    end

    H = MPO(ampo, sites)

    return sites, H
end

function main()
    args = get_args("tv")
    @show args["out_filename"]

    sites, H = get_H(args)
    psi = get_init_psi_fermion(args, sites)
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
