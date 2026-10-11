# Compile into a temporary directory for the whole test run, so the tests never
# read or write the user's cache of compiled models. Both calls are scoped to
# teardown_env(): with the default scope they would be undone as soon as this
# file has been sourced.
withr::local_envvar(
  BAYESEFRON_CACHE_ROOT = withr::local_tempdir(.local_envir = teardown_env()),
  .local_envir = teardown_env()
)

# The message that announces the new default grid is shown once in a
# session. It is marked as shown here, so that it does not appear in the
# output of whichever test happens to fit first; test-fit-settings.R clears
# the mark where it tests the message.
local({
  shown <- getFromNamespace(".bayesEfron_msgs_emitted", "bayesEfron")
  shown[["default_grid_0.4"]] <- TRUE
  # The same for the warning that `backend = "base"` of plot() is deprecated;
  # test-plot-base.R clears the mark where it tests the warning.
  shown[["plot_backend_base"]] <- TRUE
})
