"""
Tests for the second difficulty axis: structure_score / difficulty_profile.

interestingness_score ranks "deep narrow wells in a flat sea" and is close to
blind to "dense multimodality on a curved surface" — over the 73-entry corpus
`corr(curvature_score, interestingness) = -0.534`, and curvature is negatively
correlated with every other term in that composite (-0.15 to -0.66), so no
reweighting of an all-positive linear form can fix it. structure_score scores
the second axis directly.

These tests pin the properties that make the second axis worth having:
it must RISE with curvature (the thing the old score inverted), FALL with
flatness (which the old score rewarded), and stay in [0,1].
"""

using Test
using DynamicObjectives
using DynamicObjectives: structure_score, difficulty_profile, interestingness_score,
                         CURVATURE_SCALE, GridScoreResult

# Build a GridScoreResult varying only the fields the structure axis reads.
function _gs(; n_local_minima = 1, curvature_score = 0.0, plateau_fraction = 0.0,
             n_deceptive = 0, dynamic_range = 0.0, basin_fraction = 1.0,
             basin_depth_ratio = 0.0, conditioning = 0.0, ruggedness = 0.0)
    return GridScoreResult(
        zeros(5, 5),
        [collect(0.0:4.0), collect(0.0:4.0)],
        [5, 5],
        n_local_minima,
        [CartesianIndex(3, 3)],
        [0.0],
        [[2.0, 2.0]],
        dynamic_range,
        basin_fraction,
        n_deceptive,
        Vector{Float64}[],
        25, 25, 10.0,
        zeros(2),
        curvature_score,
        basin_depth_ratio,
        conditioning,
        ruggedness,
        plateau_fraction,
    )
end

@testset "structure_score / difficulty_profile" begin

    @testset "stays in [0, 1] across the corner cases" begin
        for gs in (
            _gs(),
            _gs(n_local_minima = 200, curvature_score = 1.0, n_deceptive = 99),
            _gs(plateau_fraction = 1.0),
            _gs(curvature_score = 1.0, plateau_fraction = 0.0, n_local_minima = 50,
                n_deceptive = 20),
        )
            s = structure_score(gs)
            @test 0.0 <= s <= 1.0
        end
    end

    @testset "RISES with curvature — the inversion this axis exists to fix" begin
        flat = structure_score(_gs(curvature_score = 0.0))
        mid = structure_score(_gs(curvature_score = 0.1))
        curved = structure_score(_gs(curvature_score = 1.0))
        @test flat < mid < curved
        # Curvature carries 0.40, and CURVATURE_SCALE rescales so that a raw
        # curvature_score of 0.02/0.1 = 0.2 already saturates the term. Anything
        # at or above that ceiling contributes the full 0.40.
        @test curved - flat ≈ 0.40 atol = 1e-9
        @test structure_score(_gs(curvature_score = 0.2)) ≈ curved atol = 1e-9
    end

    @testset "FALLS with flatness — which interestingness_score rewards" begin
        sharp = structure_score(_gs(plateau_fraction = 0.0))
        flat = structure_score(_gs(plateau_fraction = 0.7))
        @test sharp > flat
        @test sharp - flat ≈ 0.20 atol = 1e-9
        # Beyond the 0.7 saturation point it cannot fall further.
        @test structure_score(_gs(plateau_fraction = 1.0)) ≈ flat atol = 1e-9

        # The contrast that motivates the score, and the reason it is a
        # DISTRIBUTIONAL defect rather than a sign error in the formula.
        #
        # At a corpus-REALISTIC curvature (median curvature_score ≈ 0.089), the
        # curvature term contributes 0.13 * 0.089 ≈ 0.012 while plateau
        # contributes its full 0.08 — so interestingness_score prefers the FLAT
        # landscape, and structure_score prefers the CURVED one. Opposite order,
        # which is exactly the -0.534 correlation seen on the real corpus.
        flat_real = _gs(plateau_fraction = 0.9, curvature_score = 0.0, n_local_minima = 20)
        curved_real =
            _gs(plateau_fraction = 0.0, curvature_score = 0.089, n_local_minima = 20)
        @test interestingness_score(flat_real) > interestingness_score(curved_real)
        @test structure_score(flat_real) < structure_score(curved_real)

        # At a SYNTHETIC curvature of 1.0 the ordering flips back, because then
        # curvature finally out-earns plateau (0.13 > 0.08). Asserted so nobody
        # "fixes" this by reading the formula alone: the arithmetic is fine, it
        # is the realized range of curvature_score that makes it inert.
        curved_max =
            _gs(plateau_fraction = 0.0, curvature_score = 1.0, n_local_minima = 20)
        @test interestingness_score(curved_max) > interestingness_score(flat_real)
        @test structure_score(curved_max) > structure_score(flat_real)
    end

    @testset "rises with multimodality and deceptiveness" begin
        @test structure_score(_gs(n_local_minima = 15)) >
              structure_score(_gs(n_local_minima = 1))
        @test structure_score(_gs(n_deceptive = 5)) > structure_score(_gs(n_deceptive = 0))
        # Multimodality saturates at 15 minima, same as the other composite.
        @test structure_score(_gs(n_local_minima = 15)) ≈
              structure_score(_gs(n_local_minima = 40)) atol = 1e-9
    end

    @testset "weights sum to 1, so a maximal landscape scores exactly 1" begin
        best = _gs(curvature_score = 1.0, n_local_minima = 15, plateau_fraction = 0.0,
                   n_deceptive = 5)
        @test structure_score(best) ≈ 1.0 atol = 1e-9
        worst = _gs(curvature_score = 0.0, n_local_minima = 1, plateau_fraction = 1.0,
                    n_deceptive = 0)
        @test structure_score(worst) ≈ 0.0 atol = 1e-9
    end

    @testset "difficulty_profile labels the four quadrants" begin
        # Deep wells in a flat sea: high resolution, low structure.
        deep = _gs(n_local_minima = 30, plateau_fraction = 1.0, dynamic_range = 6.0,
                   basin_fraction = 0.0, basin_depth_ratio = 1.0, conditioning = 4.0,
                   ruggedness = 1.0, n_deceptive = 5, curvature_score = 0.0)
        p = difficulty_profile(deep)
        @test p.resolution >= 0.5
        @test p.structure < 0.5
        @test p.mode == :deep_well

        # Curved multimodal: the goodwin4d shape.
        curved = _gs(n_local_minima = 40, curvature_score = 1.0, plateau_fraction = 0.05,
                     n_deceptive = 3)
        p2 = difficulty_profile(curved)
        @test p2.structure >= 0.5
        @test p2.mode in (:curved, :both)

        # Flat, unimodal, no range: easy on both axes.
        p3 = difficulty_profile(_gs())
        @test p3.mode == :easy

        # The profile must agree with the standalone functions.
        @test p.resolution == interestingness_score(deep)
        @test p.structure == structure_score(deep)
    end

    @testset "CURVATURE_SCALE is the documented corpus calibration" begin
        # 0.02 is the corpus 75th percentile of min(directional_variation).
        # The old /0.1 normalization left 72 of 73 entries unable to saturate;
        # if this constant drifts, the axis silently loses its range again.
        @test CURVATURE_SCALE ≈ 0.02
        @test 0.0 < CURVATURE_SCALE < 0.1
    end
end
