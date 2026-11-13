# CanopyExchange maintenance

Keep this package dependency-free and independent of leaf physiology, radiation,
thermodynamics, host fields, turbulence and time integration. It owns area/volume
measures and conservative transfer arithmetic, not the physical exchange laws.
Preserve Float32/Float64 and explicit absent oxygen. Parameters and tolerances
belong in TOML fixtures. Major steps need independent code and physics reviews.
Keep runtime Project.toml dependency-free; do not add private dependencies.
