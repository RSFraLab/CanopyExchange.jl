# CanopyExchange maintenance

Keep this package free of runtime dependencies and independent of leaf physiology,
radiation, thermodynamics, host fields, turbulence and time integration. It owns
area/volume measures and conservative transfer arithmetic, not the physical exchange
laws. Preserve Float32/Float64 and the tracer rule: a tracer that is not computed is
not a key of the tracer `NamedTuple`, never a stored zero; tracers are molar and
passive. Parameters and tolerances belong in TOML fixtures. Never regenerate
`test/reference/exchange_v0_1_1.toml` from current code. Major steps need independent
code and physics reviews. Keep runtime `[deps]` empty; ForwardDiff is only a weak
dependency for the scalar-hook extension; do not add private dependencies.
