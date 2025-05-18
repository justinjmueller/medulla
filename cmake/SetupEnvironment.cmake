# cmake/SetupEnvironment.cmake

# Set paths to external dependencies based on environment variables

set(SBNANA_BASE "$ENV{MRB_INSTALL}/sbnana/$ENV{SBNANA_VERSION}")
set(SBNANAOBJ_BASE "$ENV{MRB_INSTALL}/sbnanaobj/$ENV{SBNANAOBJ_VERSION}")

set(SRPROXY_INC "$ENV{SRPROXY_INC}")
set(OSCLIB_INC "$ENV{OSCLIB_INC}")
set(OSCLIB_LIB "$ENV{OSCLIB_LIB}")
set(EIGEN_INC "$ENV{EIGEN_INC}")

set(SBNANA_INC "${SBNANA_BASE}/include" CACHE INTERNAL "")
set(SBNANA_LIB "${SBNANA_BASE}/slf7.x86_64.e26.prof/lib" CACHE INTERNAL "")
set(SBNANAOBJ_INC "${SBNANAOBJ_BASE}/include" CACHE INTERNAL "")
set(SBNANAOBJ_LIB "${SBNANAOBJ_BASE}/slf7.x86_64.e26.prof/lib" CACHE INTERNAL "")