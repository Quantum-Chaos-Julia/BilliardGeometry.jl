# Package-level quality checks with Aqua, one testset per check, on the module BilliardGeometry:
#   src/BilliardGeometry.jl   (export list and module structure) and Project.toml
# Checks: undefined_exports, unbound_args, ambiguities (recursive=false, only methods owned by BilliardGeometry),
# stale_deps, deps_compat, project_extras, persistent_tasks, piracies. undocumented_names is not run.
# Plus two export-hygiene checks: `angle` must not be exported (it would shadow Base.angle) and the undefined name
# `LimaconSegment` must not be exported. No eigen solves, no pinned wavefunction values.

using Aqua

@testset "Aqua undefined_exports" begin
    # G10: LimaconSegment is exported but never defined
    Aqua.test_undefined_exports(BilliardGeometry)
end

@testset "Aqua unbound_args" begin
    Aqua.test_unbound_args(BilliardGeometry)
end

@testset "Aqua ambiguities" begin
    Aqua.test_ambiguities(BilliardGeometry; recursive=false)
end

@testset "Aqua stale_deps" begin
    Aqua.test_stale_deps(BilliardGeometry)
end

@testset "Aqua deps_compat" begin
    Aqua.test_deps_compat(BilliardGeometry)
end

@testset "Aqua project_extras" begin
    Aqua.test_project_extras(BilliardGeometry)
end

@testset "Aqua persistent_tasks" begin
    Aqua.test_persistent_tasks(BilliardGeometry)
end

@testset "Aqua piracies" begin
    Aqua.test_piracies(BilliardGeometry)
end

@testset "export hygiene (angle, LimaconSegment)" begin
    # G10: `export angle` shadows Base.angle, so `angle(a, b)` after `using BilliardGeometry` is ambiguous in Main
    @test !(:angle in names(BilliardGeometry))
    # G10: LimaconSegment is not defined anywhere in src/, so it must not be exported
    @test !(:LimaconSegment in names(BilliardGeometry))
end
