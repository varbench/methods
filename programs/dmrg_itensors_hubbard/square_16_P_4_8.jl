using ITensors, LinearAlgebra

function snake_order(nx,ny)
    order = zeros(Int64,nx,ny)
    count = 1
    for i in 1:nx
        for j in 1:ny
            order[i,j] = count
            count += 1
        end
        if i%2 == 0
            order[i,:] =  reverse(order[i,:],dims=1)
        end
    end
    return order
end

function chain_edges(N,periodic)
    edges = []
    for i in 1:N-1
        edge = [i, i+1]
        insert!(edges,1,edge)
    end
    if periodic
        edge = [1,N]
        insert!(edges,1,edge)
    end
    return edges
end

function square_edges(order,xperiodic,yperiodic)
    nx = size(order)[1]
    ny = size(order)[2]
    edges = []
    for i in 1:nx, j in 1:ny
        for istep in -1:1, jstep in -1:1
            if (istep != 0 || jstep != 0) && (istep == 0 || jstep == 0)
                try
                    edge1 = [order[i+istep,j+jstep],order[i,j]]
                    edge2 = [order[i,j],order[i+istep,j+jstep]]
                    !(edge1 in edges) && insert!(edges,1,edge2)
                catch
                end
            end
        end
    end
    if xperiodic
        for i in 1:ny
            edge = [order[1,i],order[nx,i]]
            insert!(edges,1,edge)
        end
    end
    if yperiodic
        for j in 1:nx
            edge = [order[j,1],order[j,ny]]
            insert!(edges,1,edge)
        end
    end
    return edges
end

function chain_hubbard_MPO(edges,sites,N,t,U)
    ampo = OpSum()
    for edge in edges
        ampo += (-t,"Cdagup",edge[1],"Cup",edge[2])
        ampo += (-t,"Cdagup",edge[2],"Cup",edge[1])
        ampo += (-t,"Cdagdn",edge[1],"Cdn",edge[2])
        ampo += (-t,"Cdagdn",edge[2],"Cdn",edge[1])
    end
    for i in 1:N
        ampo += U,"Nupdn",i
    end
    H = MPO(ampo,sites)
    return H
end

function lattice_hubbard_MPO(edges,sites,order,t,U)
    nx = size(order)[1]
    ny = size(order)[2]
    N = nx*ny
    ampo = OpSum()

    for edge in edges
        ampo += (-t,"Cdagup",edge[1],"Cup",edge[2])
        ampo += (-t,"Cdagup",edge[2],"Cup",edge[1])
        ampo += (-t,"Cdagdn",edge[1],"Cdn",edge[2])
        ampo += (-t,"Cdagdn",edge[2],"Cdn",edge[1])
    end

    for i in 1:N
        ampo += U,"Nupdn",i
    end

    H = MPO(ampo,sites)
    return H
end

let
    U = 8
    nx = 4
    ny = 4
    N = nx*ny
    yperiodic = true
    xperiodic = true

    site_type = "Electron"
    sites = siteinds(site_type,N;conserve_qns=true)

    spins = [isodd(i) ? "Up" : "Dn" for i in 1:N]

    for i in [i for i in 5:12]
        spins[i] = "Emp"
    end
    @show count(i->i=="Up",spins)

    ψ = productMPS(sites,spins)

    order = snake_order(nx,ny)
    edges = square_edges(order,xperiodic,yperiodic)
    H = lattice_hubbard_MPO(edges,sites,order,1,U)

    dims = [4,8,16,32,64,128,256,500,800,1200,1700,2400,3200]
    for dim in dims
        sweeps = Sweeps(2)
        setmaxdim!(sweeps,dim)

        E,ψ = dmrg(H,ψ,sweeps)
        @show E

        var = inner(H,ψ,H,ψ) - E^2
        @show var
    end
end
