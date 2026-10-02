# Compile into a temporary directory for the whole test run, so the tests never
# read or write the user's cache of compiled models. Both calls are scoped to
# teardown_env(): with the default scope they would be undone as soon as this
# file has been sourced.
withr::local_envvar(
  BAYESEFRON_CACHE_ROOT = withr::local_tempdir(.local_envir = teardown_env()),
  .local_envir = teardown_env()
)
