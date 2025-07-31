# cmake/CheckEnvironment.cmake

macro(require_env_var VAR_NAME)
    if(NOT DEFINED ENV{${VAR_NAME}} OR "$ENV{${VAR_NAME}}" STREQUAL "")
        message(FATAL_ERROR "Required environment variable '${VAR_NAME}' is not set.")
    endif()
endmacro()

# Required for path setup
require_env_var(SBNANA_VERSION)
require_env_var(SBNANAOBJ_VERSION)

# Optional: warn if unset but don't fail
foreach(opt_var SRPROXY_INC OSCLIB_INC OSCLIB_LIB EIGEN_INC)
    if(NOT DEFINED ENV{${opt_var}} OR "$ENV{${opt_var}}" STREQUAL "")
        message(WARNING "Optional environment variable '${opt_var}' is not set.")
    endif()
endforeach()