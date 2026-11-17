# Frozen reference record

`exchange_v0_1_1.toml` holds IEEE-754 bit strings (`bitstring`) of CanopyExchange
outputs written by `test/write_reference_record.jl` from **v0.1.1, commit `c35ad8d`**,
with Julia 1.12.6, before any refactored code existed. It is compared bit for bit by
`test/test_reference_record.jl`. Never regenerate it from later code; write comparison
output to another path.

Inputs come from `test/parameters/conservative_exchange.toml`. There are four
profiles: `Float32` and `Float64`, each without and with O₂. Each entry is a list (one
element per layer, face or stream) of component lists in the order of the v0.1.x
records:

- exchange records: `(energy, water, carbon[, oxygen])`, now
  `(energy, water, tracers.co2[, tracers.o2])`;
- moist-air records: `(rho_total, rho_water, rho_energy, carbon[, oxygen])`.

| Key | Quantity | Units of (energy, water, carbon, oxygen) |
| --- | --- | --- |
| `leaf_flux`, `shifted_leaf_flux` | leaf fluxes, before and after the datum shift | W m⁻² leaf, kg m⁻² leaf s⁻¹, mol m⁻² leaf s⁻¹ |
| `dry_stem_flux`, `wet_stem_flux`, `shifted_stem_flux` | stem fluxes | per m² stem |
| `shifted_face_flux`, `steady_top_flux`, `accumulated_steady_fluxes`, `surface_static_energy_flux` | face fluxes, upward | W m⁻², kg m⁻² s⁻¹, mol m⁻² s⁻¹ (ground) |
| `leaf_air_source`, `leaf_donor_source`, `dry_plant_*_source`, `wet_plant_*_source`, `shifted_volume_source` | volume sources | W m⁻³, kg m⁻³ s⁻¹, mol m⁻³ s⁻¹ (bulk) |
| `layer_tendency`, `plant_layer_tendency` | bulk-volume tendencies | W m⁻³, kg m⁻³ s⁻¹, mol m⁻³ s⁻¹ (bulk) |
| `moist_air_unit_gas_fraction`, `moist_air_gas_fraction` | gas-volume tendencies at ε = 1 and ε = 0.96 | kg m⁻³ s⁻¹ (total, water), W m⁻³, mol m⁻³ s⁻¹ (gas) |
