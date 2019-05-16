using Printf

struct Atoms
    names::Array{String, 1}
    natom::Int64
    atom_types::Array{Int64, 1}
    acell::Array{Vector{Float64}, 1}
    cart_coors::Array{Vector{Float64}, 1}
    red_coors::Array{Vector{Float64}, 1}
end

function red2cart(red_coors::Array{Vector{Float64}, 1}, acell::Array{Vector{Float64}, 1})
    #= red coordinates -> cart coordinates =#
    natom = length(red_coors)
    acell_matrix = vcat(map(x -> x', acell)...)
    red_coors_matrix = vcat(map(x -> x', red_coors)...)
    cart_coors = red_coors_matrix * acell_matrix
    cart_coors = [cart_coors[index, :] for index in 1:natom]
    return cart_coors
end

function cart2red(cart_coors::Array{Vector{Float64}, 1}, acell::Array{Vector{Float64}, 1})
    #= cart coordinates -> red coordinates =#
    natom = length(cart_coors)
    acell_matrix = vcat(map(x -> x', acell)...)
    cart_coors_matrix = vcat(map(x -> x', cart_coors)...)
    cart_coors = cart_coors_matrix * inv(acell_matrix)
    cart_coors = [cart_coors[index, :] for index in 1:natom]
    return cart_coors
end


function poscar2Atoms(filename)
    lines = open(filename, "r") do f; readlines(f) end
    scale = parse(Float64, lines[2])
    acell = map(x -> parse.(Float64, x), split.(lines[3:5]))
    names = split(lines[6])
    n_atoms = parse.(Int64, split(lines[7]))
    natom = sum(n_atoms)
    atom_types = vcat([repeat([index], n) for (index, n) in enumerate(n_atoms)]...)
    is_selective::Bool = lines[8][1] in "sS"
    coordinate_type::Char = lines[8+is_selective][1]
    coordinates = map(x -> parse.(Float64, x[1:3]), split.(lines[9+is_selective:8+natom+is_selective]))
    if coordinate_type in "dD"
        cart_coors = red2cart(coordinates, scale*acell)
        red_coors = coordinates
    elseif coordinate_type in "cC"
        cart_coors = coordinates
        red_coors = cart2red(coordinates, scale*acell)
    end
    @assert length(coordinates) == natom
    return Atoms(names, natom, atom_types, scale*acell, cart_coors, red_coors)
end

function generate_salmon(atoms::Atoms, filename)
    #=     
    generate salmon.input file
    =#
    names = atoms.names
    pseudo_files = map(x -> "./psps/$(x)", filter(x -> occursin(".fhi", x), readdir("./psps/")))
    open(filename, "w") do f
        # angst, eV, fs
        @printf(f, "&units\n")
        @printf(f, "  unit_system = 'A_eV_fs'\n")
        @printf(f, "/\n\n")

        # calculation
        @printf(f, "&calculation\n")
        @printf(f, "  calc_mode = 'GS'\n")
        @printf(f, "/\n\n")

        # controll
        @printf(f, "&control\n")
        @printf(f, "  sysname = 'test'\n")
        @printf(f, "/\n\n")

        # rgrid
        @printf(f, "&rgrid\n")
        @printf(f, "  dl = 0.10, 0.10, 0.10\n")
        # @printf(f, "  num_rgrid = 12, 12, 12\n")
        @printf(f, "/\n\n")

        # scf
        @printf(f, "&scf\n")
        @printf(f, "  ncg = 4\n")
        @printf(f, "  nscf = 1000\n")
        @printf(f, "  convergence = 'norm_rho_dng'\n")
        @printf(f, "  threshold_norm_rho = 1.0e-15\n")
        @printf(f, "/\n\n")

        # functional
        @printf(f, "&functional\n")
        @printf(f, "  xc = 'PZ'\n")
        @printf(f, "/\n\n")

        # analysis
        @printf(f, "&analysis\n")
        @printf(f, "  out_psi = 'y'\n")
        @printf(f, "  out_dos = 'y'\n")
        @printf(f, "  out_pdos = 'y'\n")
        @printf(f, "  out_dns = 'y'\n")
        @printf(f, "  out_elf = 'y'\n")
        @printf(f, "/\n\n")

        # pseudo
        nelec = 0 # number of electrons
        @printf(f, "&pseudo\n")
        for (index, name) in enumerate(names)
            filename_psp = pseudo_files[index]
            @printf(f, "  pseudo_file(%d) = '%s'\n", index, filename_psp)
            open(filename_psp, "r") do f_psp
                for (index_line, line) in enumerate(eachline(f_psp))
                    if index_line == 2
                        zatom, zion = parse.(Float64, split(line)[1:2])
                        nelec += zion
                        @printf(f, "  izatom(%d) = %d\n", index, zatom)
                    elseif index_line == 3
                        lmax, lloc = parse.(Int64, split(line)[3:4])
                        @printf(f, "  lmax_ps(%d) = %d\n", index, lmax)
                        @printf(f, "  lloc_ps(%d) = %d\n", index, lloc)
                        break
                    end
                end
            end
        end
        @printf(f, "/\n\n")

        # system
        @printf(f, "&system\n")
        @printf(f, "  iperiodic = 0\n")
        @printf(f, "  al = % .4e, % .4e, % .4e\n", atoms.acell[1][1], atoms.acell[2][2], atoms.acell[3][3])
        @printf(f, "  nstate = %d\n", nelec)
        @printf(f, "  nelem = %d\n", length(atoms.names))
        @printf(f, "  natom = %d\n", atoms.natom)
        @printf(f, "  nelec = %d\n", nelec)
        @printf(f, "/\n\n")


        @printf(f, "&atomic_coor\n")
        for (coordinate, atom_id) in zip(atoms.cart_coors, atoms.atom_types)
            @printf(f, "  '%s' % .4f % .4f % .4f %d\n", names[atom_id], coordinate..., atom_id)
        end
        @printf(f, "/\n\n")


        # parallel
        # @printf(f, "&parallel\n")
        # @printf(f, "  nproc_k = 1\n")
        # @printf(f, "  nproc_ob = 1\n")
        # @printf(f, "  nproc_domain = 1 1 1\n")
        # @printf(f, "  nproc_domain_s = 1 1 1\n")
        # @printf(f, "/\n")
    end
end

function main()
    filename = "POSCAR"
    atoms = poscar2Atoms(filename)
    generate_salmon(atoms, "salmon.inp")
end

main()