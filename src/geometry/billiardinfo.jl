################################################################################
# BILLIARD CONSTRUCTION SUMMARY
#
# `show` of a billiard (and `log_billiard`) describes how it is built, from the
# fundamental domain down to the boundary conditions of its curves, followed by
# the symmetry registry, the symmetry group G' and a few global geometric
# quantities. Curve parameters are not shown.
################################################################################

_sym_label(g::AbsSymmetry) = string(g)
_sym_label(g::NFoldRotation) = "NFoldRotation(2π*$(g.m)/$(g.order), $(g.sym_id))"

_bc_text(::SpecularReflection, reg) = "SpecularReflection (physical wall)"
_bc_text(::QuantumSolverIgnore, reg) = "QuantumSolverIgnore (virtual edge, no condition imposed)"
_bc_text(bc::Transparent, reg) = "Transparent (join to domain $(bc.next_id))"
_bc_text(bc::PeriodicX, reg) = "PeriodicX (linked to domain $(bc.next_id))"
function _bc_text(bc::SymmetryWall, reg)
    i = findfirst(s -> s.sym_id == bc.sym_id, reg)
    label = isnothing(i) ? "sym_id $(bc.sym_id)" : _sym_label(reg[i])
    return "SymmetryWall (cut of $label, neighbouring sector $(bc.sector_id))"
end
_bc_text(bc::AbsBoundaryCondition, reg) = string(nameof(typeof(bc)))

_hole_registry(::AbsDomainHole) = SymmetryRegistry(AbsSymmetry[])
_hole_registry(h::SymmetricHole) = h.symmetries
_hole_registry(h::SymmetrizedHoles) = _hole_registry(h.base)

function _write_domain!(io::IO, d::AbsSimpleDomain, reg, ind::Int)
    println(io, " "^ind, "domain ", d.id, " (", nameof(typeof(d)), "):")
    for c in d.boundary
        println(io, " "^(ind + 2), "curve ", c.segment_id, " (", nameof(typeof(c)), "): ", _bc_text(c.bc, reg))
    end
end

function _group_label(group)
    order = length(group) + 1
    order == 1 && return "trivial"
    nrot = 1 + count(g -> !orientation_reversing(g), group)
    return nrot == order ? "C_$nrot" : "D_$nrot"
end

_hole_center(::AbsDomainHole) = nothing
_hole_center(h::SymmetricHole) = Tuple(h.center)
_hole_center(h::SymmetrizedHoles) = _hole_center(h.base)

_symmetry_keys(reg) = Set(symmetry_key(g) for g in reg)

_fundamental_hole_domains(h::AbsDomainHole) = get_all_domains(h)
_fundamental_hole_domains(h::SymmetricHole) = get_all_domains(h.piece)
_fundamental_hole_domains(h::SymmetrizedHoles) = _fundamental_hole_domains(h.base)

_outer_subdomains_of(d::AbsMultiplyConnectedDomain) = d.subdomains
_outer_subdomains_of(d::AbsDomain) = get_all_domains(d)
_holes_of(::AbsDomain) = AbsDomainHole[]
_holes_of(d::AbsMultiplyConnectedDomain) = d.holes

function _corner_summary(billiard::AbsBilliard)
    angles = corner_angles(billiard; fundamental=false)
    isempty(angles) && return "none (smooth boundary)"
    fracs = round.(angles ./ pi; digits=4)
    counts = [(a, count(==(a), fracs)) for a in sort!(unique(fracs))]
    return "$(length(angles)): " * join(("$(c)×$(a)π" for (a, c) in counts), ", ")
end

function _write_boundary_length!(io::IO, billiard::AbsBilliard)
    map = full_boundary_components(billiard)
    lengths = [last(component_lengths(loop)) for loop in map.curves]
    if length(lengths) == 1
        println(io, "    boundary length: ", lengths[1])
        return
    end
    println(io, "    boundary length: ", sum(lengths))
    for (i, (L, hole)) in enumerate(zip(lengths, map.is_hole))
        println(io, "      ", hole ? "hole loop $(i - 1)" : "outer boundary", ": ", L)
    end
end

function Base.show(io::IO, ::MIME"text/plain", billiard::AbsBilliard)
    fd = billiard.fundamental_domain
    println(io, nameof(typeof(billiard)), ":")
    println(io, "  Fundamental domain: ", nameof(typeof(fd)))
    println(io, "  Domains and boundary conditions:")
    for d in _outer_subdomains_of(fd)
        _write_domain!(io, d, billiard.symmetries, 4)
    end
    holes = _holes_of(fd)
    if !isempty(holes)
        println(io, "  Holes:")
        for (i, h) in enumerate(holes)
            println(io, "    hole ", i, " (", nameof(typeof(h)), ", ", hole_loop_count(h), " boundary loop(s)):")
            for d in _fundamental_hole_domains(h)
                _write_domain!(io, d, _hole_registry(h), 6)
            end
        end
    end
    group = get_symmetry_group(billiard)
    println(io, "  Symmetries (", _group_label(group), "):")
    if isempty(billiard.symmetries)
        println(io, "    none")
    else
        println(io, "    registered elements:")
        for g in billiard.symmetries
            println(io, "      ", _sym_label(g))
        end
        gens = get_symmetries(billiard)
        isempty(gens) || println(io, "    generators: ", join(_sym_label.(gens), ", "))
    end
    main_keys = _symmetry_keys(billiard.symmetries)
    differing = [(i, h) for (i, h) in enumerate(holes) if !isempty(_hole_registry(h)) && _symmetry_keys(_hole_registry(h)) != main_keys]
    if !isempty(differing)
        println(io, "  Hole Symmetries:")
        for (i, h) in differing
            reg = _hole_registry(h)
            println(io, "    hole ", i, " (", _group_label(collect(reg)), ") about ", _hole_center(h), ":")
            for g in reg
                println(io, "      ", _sym_label(g))
            end
        end
    end
    placed = [(i, h) for (i, h) in enumerate(holes) if h isa SymmetrizedHoles]
    if !isempty(placed)
        println(io, "  Symmetrized Hole:")
        for (j, h) in placed
            println(io, "    hole ", j, " (", length(h.holes), " copies of the base hole, placed by):")
            for g in h.symmetries
                println(io, "      ", _sym_label(g))
            end
        end
    end
    println(io, "  Geometry:")
    println(io, "    holes (genus): ", genus(billiard))
    println(io, "    area: ", area(billiard))
    println(io, "    fundamental area: ", fundamental_area(billiard))
    _write_boundary_length!(io, billiard)
    print(io, "    corners: ", _corner_summary(billiard))
end

"""
    log_billiard(billiard::AbsBilliard) → nothing

Log with `@info` how `billiard` is constructed: its domains and the boundary
condition of each curve (`domain id → curve id → boundary condition`), its holes,
the registered symmetry elements, the symmetry group `G'` and its generators, and
the area, fundamental area, genus and corners. Curve parameters are not shown.
The same text is displayed by `show` when a billiard is returned at the REPL.
"""
function log_billiard(billiard::AbsBilliard)
    @info repr(MIME("text/plain"), billiard)
    return nothing
end
