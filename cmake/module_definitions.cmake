#
# Define variables of module groups for use by the sebsjames/oces_viewer
# build process itself, and by client projects.
#

macro(setup_module_variables_for_oces_reader oces_directory maths_directory)

  include(${maths_directory}/cmake/module_definitions.cmake)
  setup_module_variables_for_maths (${maths_directory} "") # Don't use json, so pass empty

  # We import sm.mathconst, sm.vec and sm.vvec as well as sm.mat and sm.geometry.
  set(OCES_READER_MATHS_MODULES
    ${SM_MATHCONST_MODULES}
    ${SM_VVEC_MODULES}
    ${SM_VEC_MODULES}
    ${SM_MAT_MODULES}
    ${maths_directory}/sm/centroid.cppm
    ${SM_GEOMETRY_MODULES}
  )
  list(REMOVE_DUPLICATES OCES_READER_MATHS_MODULES)

  set(OCES_READER_MODULES
    ${OCES_READER_MATHS_MODULES}
    ${oces_directory}/oces/reader.cppm
  )
  list(REMOVE_DUPLICATES OCES_READER_MODULES)

  set(OCES_HEXEYE_MODULES
    ${OCES_READER_MODULES}
    ${oces_directory}/oces/hexeye.cppm
    ${SM_VVEC_MODULES}
    ${SM_HEXGRID_HDF_MODULES}
  )
  list(REMOVE_DUPLICATES OCES_HEXEYE_MODULES)

endmacro()
