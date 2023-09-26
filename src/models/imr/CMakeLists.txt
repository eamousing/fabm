add_library(fabm_models_imr OBJECT
    imr_model_library.F90
    norwecom_num.F90
    num_model/Fortran/read_input.f90
    num_model/Fortran/globals.f90
    num_model/Fortran/spectrum.f90
    num_model/Fortran/generalists_simple.f90
    num_model/Fortran/generalists.f90
    num_model/Fortran/diatoms_simple.f90
    num_model/Fortran/diatoms.f90
    num_model/Fortran/copepods.f90
    num_model/Fortran/POM.f90
    num_model/Fortran/NUMmodel.f90
    )

# Run the fypp preprocessor
# Code from https://fypp.readthedocs.io/en/stable/fypp.html#examples
set(fppFiles norwecom_num.fpp)

foreach(infileName IN LISTS fppFiles)
    # Generate output file name
    string(REGEX REPLACE ".fpp\$" ".F90" outfileName "${infileName}")

    # Create the full path for the new file
    set(outfile "${CMAKE_CURRENT_BINARY_DIR}/${outfileName}")

    # Generate input file name
    set(infile "${CMAKE_CURRENT_SOURCE_DIR}/${infileName}")

    # Custom command to do the processing
    add_custom_command(
        OUTPUT "${outfile}"
        COMMAND fypp "${infile}" "${outfile}"
        MAIN_DEPENDENCY "${infile}"
        VERBATIM)

    # Finally add output file to a list
    set(outFiles ${outFiles} "${outfile}")
endforeach(infileName IN LISTS fppFiles)

add_dependencies(fabm_models_imr fabm_base)