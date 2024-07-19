#!/usr/bin/env julia
#
# 1D t-V model

include("common.jl")

function get_H(args)
    @assert args["L2"] == args["L"]
    @assert !args["mars"]
    @assert args["J2"] == 0
    @assert args["J22"] == 0
    @assert args["h"] == 0

    L = args["L"]
    peri = args["peri"]
    V = args["V"]

    sites = siteinds("Fermion", L; conserve_qns = true)

    ampo = OpSum()

    function add_edge!(i1, i2)
        i1 = mod1(i1, L)
        i2 = mod1(i2, L)
        ampo -= 1, "Cdag", i1, "C", i2
        ampo += 1, "C", i1, "Cdag", i2
        ampo += V, "N", i1, "N", i2
    end

    for i = 1:L-1+peri
        add_edge!(i, i + 1)
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
