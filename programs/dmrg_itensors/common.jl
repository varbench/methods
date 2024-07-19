import Pkg

Pkg.activate(".")
Pkg.instantiate()

using HDF5
using ITensors
using MKL

include("args.jl")

function get_init_psi_spin(args, sites)
    @assert args["Nf"] == 0

    dtype = args["cmpl"] ? ComplexF64 : Float64
    psi = nothing
    if args["zero_mag"]
        N_half = div(length(sites), 2)
        state = shuffle([fill("Up", N_half); fill("Dn", N_half)])
        psi = randomMPS(dtype, sites, state; linkdims = args["max_B"])
    else
        psi = randomMPS(dtype, sites; linkdims = args["max_B"])
    end
    return psi
end

function get_init_psi_fermion(args, sites)
    @assert !args["zero_mag"]

    dtype = args["cmpl"] ? ComplexF64 : Float64
    N = length(sites)
    Nf = args["Nf"]
    state = shuffle([fill("1", Nf); fill("0", N - Nf)])
    psi = randomMPS(dtype, sites, state; linkdims = args["max_B"])
    return psi
end

function get_sweeps(args)
    sweeps = Sweeps(args["max_step"])
    setmaxdim!(sweeps, args["max_B"])
    setcutoff!(sweeps, 1e-12)
    noise = 10 .^ LinRange(-3, -12, div(args["max_step"], 2))
    setnoise!(sweeps, noise..., 0)
    return sweeps
end
