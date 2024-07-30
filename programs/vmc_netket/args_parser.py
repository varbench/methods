import argparse
import os
from datetime import datetime

import numpy as np


def get_parser():
    parser = argparse.ArgumentParser(allow_abbrev=False)

    group = parser.add_argument_group("physics parameters")
    group.add_argument(
        "--ham",
        type=str,
        default="ising",
        choices=["ising", "heis", "heis_tri", "heis_kag", "j1j2", "hubb", "tv"],
        help="Hamiltonian type",
    )
    group.add_argument(
        "--boundary",
        type=str,
        default="open",
        choices=["open", "peri"],
        help="boundary conditions",
    )
    group.add_argument(
        "--sign",
        type=str,
        default="none",
        choices=["none", "mars"],
        help="sign rule",
    )
    group.add_argument(
        "--ham_dim",
        type=int,
        default=1,
        choices=[1, 2],
        help="dimension of the lattice",
    )
    group.add_argument(
        "--L",
        type=int,
        default=4,
        help="edge length of the lattice",
    )
    group.add_argument(
        "--L2",
        type=int,
        default=0,
        help="another edge length of the lattice",
    )
    group.add_argument(
        "--J2",
        type=float,
        default=0,
        help="2nd nearest neighbor interaction",
    )
    group.add_argument(
        "--U",
        type=float,
        default=0,
        help="on-site interaction",
    )
    group.add_argument(
        "--V",
        type=float,
        default=0,
        help="repulsive interaction",
    )
    group.add_argument(
        "--h",
        type=float,
        default=0,
        help="external field",
    )
    group.add_argument(
        "--zero_mag",
        action="store_true",
        help="use zero magnetization constraint",
    )
    group.add_argument(
        "--Nf",
        type=int,
        default=0,
        help="number of fermions",
    )

    group = parser.add_argument_group("network parameters")
    group.add_argument(
        "--net",
        type=str,
        default="jas",
        choices=["jas", "rbm", "gcnn", "rnn_lstm"],
        help="network type",
    )
    group.add_argument(
        "--layers",
        type=int,
        default=1,
        help="number of layers",
    )
    group.add_argument(
        "--features",
        type=int,
        default=1,
        help="number of features",
    )
    group.add_argument(
        "--dtype",
        type=str,
        default="float32",
        choices=["float32", "float64", "complex64", "complex128"],
        help="data type",
    )

    group = parser.add_argument_group("optimizer parameters")
    group.add_argument(
        "--seed",
        type=int,
        default=0,
        help="random seed, 0 for randomized",
    )
    group.add_argument(
        "--optimizer",
        type=str,
        default="sr",
        choices=["adam", "sgd", "sr"],
        help="optimizer type",
    )
    group.add_argument(
        "--split_real",
        action="store_true",
        help="split real and imaginary parts of parameters in the optimizer",
    )
    group.add_argument(
        "--batch_size",
        type=int,
        default=1024,
        help="batch size",
    )
    group.add_argument(
        "--lr",
        type=float,
        default=1e-3,
        help="learning rate",
    )
    group.add_argument(
        "--diag_shift",
        type=float,
        default=1e-2,
        help="diagonal shift of SR",
    )
    group.add_argument(
        "--max_step",
        type=int,
        default=10**4,
        help="number of training/sampling steps",
    )
    group.add_argument(
        "--grad_clip",
        type=float,
        default=0,
        help="global norm to clip gradients, 0 for disabled",
    )
    group.add_argument(
        "--chunk_size",
        type=int,
        default=1024,
        help="chunk size, 0 for disabled",
    )
    group.add_argument(
        "--estim_size",
        type=int,
        default=1024**2,
        help="batch size to estimate the Hamiltonian, 0 for matching `batch_size`",
    )

    group = parser.add_argument_group("system parameters")
    group.add_argument(
        "--show_progress",
        action="store_true",
        help="show progress",
    )
    group.add_argument(
        "--cuda",
        type=str,
        default="auto",
        help="GPU ID, empty string for disabling GPU, multi-GPU is not supported yet",
    )
    group.add_argument(
        "--run_name",
        type=str,
        default="",
        help="output subdirectory to keep repeated runs, empty string for disabled",
    )
    group.add_argument(
        "-o",
        "--out_dir",
        type=str,
        default="./out",
        help="output directory, empty string for disabled",
    )

    return parser


def get_ham_net_name(args):
    ham_name = "{ham}_{boundary}"
    if args.sign != "none":
        ham_name += "_{sign}"
    ham_name += "_{ham_dim}d_L{L}"
    if args.L2 and args.L2 != args.L:
        ham_name += ",{L2}"
    if args.J2:
        ham_name += "_J2={J2:g}"
    if args.U:
        ham_name += "_U={U:g}"
    if args.V:
        ham_name += "_V={V:g}"
    if args.h:
        ham_name += "_h={h:g}"
    if args.zero_mag:
        ham_name += "_zm"
    if args.Nf:
        ham_name += "_Nf{Nf}"
    ham_name = ham_name.format(**vars(args))

    net_name = "{net}"
    if args.net == "rbm":
        net_name += "_a{features}"
    elif args.net != "jas":
        net_name += "_l{layers}_f{features}"

    net_name += "_{optimizer}"
    if args.split_real:
        net_name += "_sp"
    if args.grad_clip:
        net_name += "_gc{grad_clip:g}"
    net_name = net_name.format(**vars(args))

    return ham_name, net_name


def post_init_args(args):
    if args.ham_dim == 1:
        assert args.L2 == 0
    else:
        if args.L2 == 0:
            args.L2 = args.L

    if args.seed == 0:
        # The seed depends on the time and the PID
        args.seed = hash((datetime.now(), os.getpid())) & (2**32 - 1)

    if args.optimizer == "sr" and args.diag_shift == 0:
        args.diag_shift = args.lr

    if args.chunk_size == 0:
        args.chunk_size = None

    if args.estim_size == 0:
        args.estim_size = args.batch_size

    args.ham_name, args.net_name = get_ham_net_name(args)

    if args.dtype in ["float32", np.float32]:
        args.dtype = np.float32
    elif args.dtype in ["float64", np.float64]:
        args.dtype = np.float64
    elif args.dtype in ["complex64", np.complex64]:
        args.dtype = np.complex64
    elif args.dtype in ["complex128", np.complex128]:
        args.dtype = np.complex128
    else:
        raise ValueError(f"Unknown dtype: {args.dtype}")

    if args.out_dir:
        args.full_out_dir = "{out_dir}/{ham_name}/{net_name}/".format(**vars(args))
        if args.run_name:
            args.full_out_dir = "{full_out_dir}{run_name}/".format(**vars(args))
        args.log_filename = args.full_out_dir + "out"
    else:
        args.full_out_dir = None
        args.log_filename = None


def set_env(args):
    if args.cuda != "auto":
        os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
        os.environ["CUDA_VISIBLE_DEVICES"] = args.cuda
    os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"

    np.random.seed(args.seed)
