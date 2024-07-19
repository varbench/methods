using ArgParse
using PyCall
using Random

function get_args(ham)
    s = ArgParseSettings()
    @add_arg_table s begin
        "--peri"
        help = "use periodic boundary conditions"
        action = :store_true
        "--mars"
        help = "use Marshall sign rule"
        action = :store_true
        "--L"
        help = "edge length of the lattice"
        arg_type = Int
        required = true
        "--L2"
        help = "another edge length of the rectangular lattice"
        arg_type = Int
        default = 0
        "--J2"
        help = "2nd nearest neighbor interaction on \\ diagonals"
        arg_type = Float64
        default = 0
        "--J22"
        help = "2nd nearest neighbor interaction on / diagonals. Set J22 = 0 for triangular lattices, and J22 = J2 for diagonal lattices"
        arg_type = Float64
        default = 0
        "--V"
        help = "repulsive interaction"
        arg_type = Float64
        default = 0
        "--h"
        help = "external field"
        arg_type = Float64
        default = 0
        "--zero_mag"
        help = "use zero magnetization constraint"
        action = :store_true
        "--Nf"
        help = "number of fermions"
        arg_type = Int
        default = 0

        "--max_B"
        help = "maximum bond dimension"
        arg_type = Int
        required = true
        "--cmpl"
        help = "use complex parameters"
        action = :store_true

        "--max_step"
        help = "number of DMRG steps"
        arg_type = Int
        default = 100
        "--seed"
        help = "random seed, 0 for randomized"
        arg_type = Int
        default = 0
    end
    args = parse_args(s)

    @assert args["L"] % 2 == 0

    out_filename = ham
    if ham != "heis_kag"
        if args["peri"]
            out_filename *= "_peri"
        else
            out_filename *= "_open"
        end
        if args["mars"]
            out_filename *= "_mars"
        end
        out_filename *= "_L{L}"
        if args["L2"] == 0
            args["L2"] = args["L"]
        else
            @assert args["L2"] % 2 == 0
            out_filename *= ",{L2}"
        end
    end
    if args["J2"] != 0
        out_filename *= "_J2={J2:g}"
    end
    if args["J22"] != 0
        out_filename *= "_J22={J22:g}"
    end
    if args["V"] != 0
        out_filename *= "_V={V:g}"
    end
    if args["h"] != 0
        out_filename *= "_J22={h:g}"
    end
    if args["zero_mag"]
        out_filename *= "_zm"
    end
    if args["Nf"] != 0
        out_filename *= "_Nf{Nf}"
    end

    out_filename *= "_B{max_B}"
    if args["cmpl"]
        out_filename *= "_cmpl"
    end

    out_filename *= ".hdf5"

    # The easiest way to do Python-like formatting is to call Python...
    py"""
def format(s, d):
    return s.format(**d)
"""
    out_filename = py"format"(out_filename, args)
    args["out_filename"] = out_filename

    if args["seed"] > 0
        Random.seed!(args["seed"])
    end

    return args
end
