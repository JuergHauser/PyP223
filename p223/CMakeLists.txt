enable_language(Fortran)
file( GLOB p223_fortran_modules *.f90 )

##add_library(p223 SHARED  ${p223_fortran_modules})

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
Python_add_library(p223 MODULE WITH_SOABI ${p223_fortran_modules})

install(TARGETS p223 DESTINATION ./pyp223)
