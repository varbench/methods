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

    periodic && insert!(edges,1,[1,N])

    return edges
end

function square_edges(order,xperiodic,yperiodic)
    nx = size(order)[1]
    ny = size(order)[2]
    edges = []
    for i in 1:nx, j in 1:ny
        for istep in -1:1, jstep in -1:1
            bool1 = (istep != 0 || jstep != 0)
            bool2 = (istep == 0 || jstep == 0)
            if bool1 && bool2
                try
                    edge = [order[i+istep,j+jstep],order[i,j]]
                    !(edge in edges) && push!(edges,reverse(edge))
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

function chain_heisenberg_MPO(edges,sites)
    ampo = OpSum()
    for edge in edges
        ampo += 2,"S+",edge[1],"S-",edge[2]
        ampo += 2,"S-",edge[1],"S+",edge[2]
        ampo += 4,"Sz",edge[1],"Sz",edge[2]
    end
    H = MPO(ampo,sites)
    return H
end

function lattice_heisenberg_MPO(NN,sites,order)
    nx = size(order)[1]
    ny = size(order)[2]
    ampo = OpSum()

    for edge in NN
        ampo += 2.0,"S+",edge[1],"S-",edge[2]
        ampo += 2.0,"S-",edge[1],"S+",edge[2]
        ampo += 4.0,"Sz",edge[1],"Sz",edge[2]
    end

    H = MPO(ampo,sites)
    return H
end

function main_function()

    #------------------------
    # HEISENBERG MODEL TESTS
    #------------------------

    #------------------------
    #CHAIN TEST
    #------------------------

    N = 200
    periodic = false

    site_type = "S=1/2"
    sites = siteinds(site_type,N;conserve_qns=true)
    ψ = productMPS(sites,n -> isodd(n) ? "Up" : "Dn")

    edges = chain_edges(N,periodic)
    H = chain_heisenberg_MPO(edges,sites)
    dims = [4,8,16,32,64,128,256,500,1000]

    for dim in dims
        sweeps = Sweeps(2)
        setmaxdim!(sweeps,dim)
        E,ψ = dmrg(H,ψ,sweeps)
        @show E

        var = inner(H,ψ,H,ψ) - E^2
        @show var
    end
end

main_function()
