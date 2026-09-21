# FastScape T coupling

These examples evolve a two-dimensional FastScape surface raster coupled to a
two-dimensional ASPECT section. The second raster direction represents the
finite transverse width of the landscape rather than a one-dimensional
transect.

Marine sediment input can be much larger than the accommodation of an
individual coastal cell during one landscape step. Set
`Limit marine deposition to available accommodation = true` to cap direct
coastal deposition at sea level plus
`Maximum marine deposition above sea level`. The solid sediment that does not
fit is conserved as exported sediment and is reported separately in
`sediment_budget.csv` as
`accommodation_limited_export_m3_per_year`. This option prevents a coastal
base-level cell from becoming an artificial sediment tower before marine
diffusion can redistribute the deposit. The default is `false` to preserve the
behavior of existing parameter files.
