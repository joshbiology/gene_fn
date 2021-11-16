# Lib

ProjectTemplate folder for any files that provide useful functionality for your work, but do not constitute a statistical analysis per se. 

`lib/helpers.R` script organizes any functions you use in your project that aren't quite general enough to belong in a package.

If you have project specific configuration that you'd like to store in the config object, you can specify that in `lib/globals.R`.  This is the first file loaded from `lib`, so any functions in `lib`, `munge` or `src` can reference this configuration by simply using the `config$my_config_var` form.
