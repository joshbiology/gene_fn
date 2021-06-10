# Add project specific configuration that can be overridden from load.project()
add.config(
  seed = 123,
  apply.override = TRUE
)

set.seed(config$seed)

#00-Settings
ggplot2::theme_set(theme_minimal())
