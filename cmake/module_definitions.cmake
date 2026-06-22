#
# Define variables of module groups for use by the sebsjames/oces_viewer
# build process itself, and by client projects.
#

# Sets up variables the define the sebsjames/maths modules used
# building oces/reader.cppm. Pass in the path to the sebsjames/maths root.
macro(setup_module_variables_for_oces_reader_maths base_directory)
  set(OCES_READER_MATHS_MODULES
    ${base_directory}/sm/mathconst.cppm;
    ${base_directory}/sm/constexpr_math.cppm;
    ${base_directory}/sm/trait_tests.cppm;
    ${base_directory}/sm/interval.cppm;
    ${base_directory}/sm/polysolve.cppm;
    ${base_directory}/sm/bessel_i0.cppm;
    ${base_directory}/sm/random.cppm;
    ${base_directory}/sm/vec.cppm;
    ${base_directory}/sm/quaternion.cppm;
    ${base_directory}/sm/mat.cppm;
    ${base_directory}/sm/algo.cppm;
    ${base_directory}/sm/vvec.cppm;
    ${base_directory}/sm/geometry.cppm;
  )
endmacro()

# Just one module at the moment
macro(setup_module_variables_for_oces_reader base_directory)
  set(OCES_READER_MODULES
    ${base_directory}/oces/reader.cppm
  )
endmacro()
