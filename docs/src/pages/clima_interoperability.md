# CliMA Interoperability

`CanopyExchange` is designed to compose directly with CliMA packages
(`Thermodynamics.jl`, `SurfaceFluxes.jl`, `ClimaLand.jl`, `ClimaAtmos.jl`, and
`ClimaCore.jl`) with no runtime dependencies in its own `Project.toml` (ForwardDiff is an optional
weak dependency that only loads the scalar-hook extension).

## 1. `SurfaceFluxes.jl` & `Thermodynamics.jl` Static-Energy Datum

In `SurfaceFluxes.jl`, surface sensible and latent heat fluxes (`sfc.shf`, `sfc.lhf`)
and kinematic evaporation (`sfc.evaporation`, $E$ [$\mathrm{kg\,m^{-2}\,s^{-1}}$])
are defined relative to `Thermodynamics.jl`'s reference temperature $T_0$ and
surface geopotential $\Phi_{\mathrm{sfc}} = g z_{\mathrm{sfc}}$:

```math
\begin{aligned}
\mathrm{SHF} &= -\rho_{\mathrm{sfc}} g_h ( \mathrm{DSE}_{\mathrm{int}} - \mathrm{DSE}_{\mathrm{sfc}} ) + \mathrm{VSE}_{\mathrm{sfc}} \, E, \\
\mathrm{LHF} &= \mathrm{LH}_{v0} \, E,
\end{aligned}
```

where $\mathrm{VSE}_{\mathrm{sfc}} = c_{pv}(T_{\mathrm{sfc}} - T_0) + \Phi_{\mathrm{sfc}}$
is `Thermodynamics.vapor_static_energy(thermo_params, T_sfc, Φ_sfc)`. Summing both
components yields the exact moist static-energy boundary flux:

```math
\mathrm{SHF} + \mathrm{LHF}
= -\rho_{\mathrm{sfc}} g_h \Delta\mathrm{DSE}
+ \bigl(h_v(T_{\mathrm{sfc}}) + \Phi_{\mathrm{sfc}}\bigr) E,
```

because $\mathrm{LH}_{v0} + c_{pv}(T_{\mathrm{sfc}} - T_0) = h_v(T_{\mathrm{sfc}})$
(`Thermodynamics.enthalpy_vapor(thermo_params, T_sfc)`). Pass `sfc.shf`, `sfc.lhf`,
and `sfc.evaporation` directly to [`surface_static_energy_flux`](@ref):

```julia
# With sfc = SurfaceFluxes.surface_conditions(sf_params, sf_inputs):
soil_flux = surface_static_energy_flux(sfc.shf, sfc.lhf, sfc.evaporation, (; co2 = soil_co2_flux))
```

For leaf and stem layers at height $z_k$ ($\Phi_k = g z_k$), pass
`Thermodynamics.enthalpy_vapor(thermo_params, T_leaf)` for a pure enthalpy budget
($\Phi = 0$), or `Thermodynamics.enthalpy_vapor(thermo_params, T_leaf) + Φ_k` when
coupling to a total moist energy budget ($\rho e_{\mathrm{tot}}$) that includes
gravitational potential energy:

```julia
hv_k = Thermodynamics.enthalpy_vapor(thermo_params, T_leaf) + grav * z_k
lf_k = sensible_water_flux(H_leaf, E_leaf, hv_k, (; co2 = C_leaf))
```

## 2. Reference-Enthalpy & Geopotential Shifts (`shift_water_energy_datum`)

Eddy-covariance flux products and standalone leaf models often express energy fluxes
relative to liquid water at $T_{\mathrm{ref}} = 273.15\,\mathrm{K}$ without
geopotential, whereas atmospheric columns use `Thermodynamics.jl`'s reference internal
energy plus local geopotential $\Phi(z)$.

Use [`shift_water_energy_datum`](@ref) with $\Delta e_w = e_{w,\mathrm{target}} - e_{w,\mathrm{source}}$
[$\mathrm{J\,kg^{-1}\,H_2O}$] on **individual shared `FaceFlux`, `LeafFlux`,
`StemFlux`, or `VolumeSource` streams prior to [`layer_tendency`](@ref)**:

```jldoctest
julia> using CanopyExchange

julia> faces = [face_flux(120.0, 4.0e-5, (; co2 = 0.0)), face_flux(150.0, 5.0e-5, (; co2 = 0.0))];

julia> z_faces = [0.0, 2.0]; grav = 9.80665;

julia> faces_geo = shift_water_energy_datum.(faces, grav .* z_faces);

julia> faces_geo[2].rates.energy - faces[2].rates.energy ≈ grav * 2.0 * 5.0e-5
true
```

Applying the shift to shared `faces` before taking finite-volume differences guarantees
that interior face fluxes cancel identically across layer interfaces; shifting a
layer's net divergence *after* `layer_tendency` would drop the vertical gradient of
$\Phi(z)$ across the cell.

## 3. Diagnostic Steady-Air Columns (`ClimaLand` Bridge)

Big-leaf and diagnostic multi-layer land models neglect canopy-air storage
($\partial_t \psi_{\mathrm{air}} = 0$) so the upward canopy-top flux equals soil
boundary flux plus vertically integrated plant sources. [`steady_top_flux`](@ref),
[`accumulate_steady_fluxes!`](@ref), and [`accumulate_steady_fluxes`](@ref) bridge
diagnostic and prognostic formulations directly:

```jldoctest
julia> using CanopyExchange

julia> geoms = [PlantLayerGeometry(1.5, 0.0, 0.2), PlantLayerGeometry(2.5, 1.4, 0.15)];

julia> lfs = [sensible_water_flux(NaN, NaN, 2.5e6, (; co2 = NaN)), sensible_water_flux(90.0, 3.5e-5, 2.5e6, (; co2 = -1.2e-5))];

julia> sfs = [sensible_stem_flux(25.0, (; co2 = 0.0)), sensible_stem_flux(15.0, (; co2 = 0.0))];

julia> air_srcs = getproperty.(plant_air_sources.(lfs, sfs, geoms), :air);

julia> soil = surface_static_energy_flux(40.0, 80.0, 3.2e-5, (; co2 = 2.0e-6));

julia> faces = accumulate_steady_fluxes(soil, air_srcs, geoms);

julia> maximum(k -> abs(layer_tendency(faces[k], faces[k + 1], air_srcs[k], geoms[k]).rates.energy), eachindex(geoms)) <= 1.0e-13
true
```

## 4. `ClimaAtmos` & `ClimaCore.Fields.Field` Porous-Canopy Tendencies

Prognostic atmospheric columns with time-invariant air porosity
$\varepsilon_{\mathrm{gas}} \in (0, 1]$ ($\partial_t \varepsilon_{\mathrm{gas}} = 0$)
and conserved dry-air mass ($S_{\rho_d} = 0$) evolve conserved densities per unit
**gas volume**. [`moist_air_tendencies`](@ref) divides bulk-volume tendencies
fieldwise by `gas_fraction`, returning [`MoistAirTendencies`](@ref) with
`rho_total == rho_water == water / ε_gas`, `rho_energy == energy / ε_gas`, and every
tracer divided by `ε_gas`. Tracers stay nested in `tracers`, so a tracer name cannot
overwrite a host field.

Every `CanopyExchange` record is a scalar in broadcasts and every kernel takes
positional arguments only, without heap allocations, so column pipelines broadcast
over `ClimaCore.Fields.Field` and `CuArray`. GPU field broadcasts reject a type or a
keyword argument in the function position, so construct records with the function
forms (`layer_geometry`, `leaf_flux`, `stem_flux`, `face_flux`) and the `Tracers`
callable (`test/run_fields.jl` runs this on CPU and CUDA fields):

```julia
# On a ClimaCore CenterFiniteDifferenceSpace:
T       = Tracers(:co2)
air(p)  = p.air
geom_f  = layer_geometry.(dz_f, lad_f, sad_f)
lf_f    = sensible_water_flux.(sh_leaf_f, evap_leaf_f, hv_f, T.(co2_leaf_f))
sf_f    = sensible_stem_flux.(sh_stem_f, T.(co2_stem_f))
air_f   = air.(plant_air_sources.(lf_f, sf_f, geom_f))
tend_f  = layer_tendency.(bot_face_f, top_face_f, air_f, geom_f)
moist_f = moist_air_tendencies.(tend_f, eps_gas_f)
```
