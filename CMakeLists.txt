###################################################################
# Package: ORStudy 
###################################################################

# declare the package name
atlas_subdir( ORStudy )

# extra deps based on the build environment:
set( extra_deps )
set( extra_libs )

# declare the ORStudy package dependencies
atlas_depends_on_subdirs(
    PUBLIC
    SusyNtuple
    ${extra_deps}
)

#external dependencies
find_package( ROOT COMPONENTS Gpad Graf Graf3d Core Tree MathCore Hist RIO )

atlas_add_root_dictionary( ORStudyLib ORStudyCintDict
        ROOT_HEADERS ORStudy/*.h Root/LinkDef.h
        EXTERNAL_PACKAGES ROOT
)

# libraries
atlas_add_library( ORStudyLib
    ORStudy/*.h Root/*.cxx ${ORStudyCintDict}
    PUBLIC_HEADERS ORStudy
    PRIVATE_INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
    LINK_LIBRARIES SusyNtupleLib
    PRIVATE_LINK_LIBRARIES ${ROOT_LIBRARIES}
)

# executable(s) in the package
set( extra_libs )

# define a function to call add_executable macro on our executables
function( SusyNtExec filename)
    set(execname)
    get_filename_component(execname ${filename} NAME_WE)
    atlas_add_executable( ${execname} "util/${execname}.cxx"
        INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
        LINK_LIBRARIES ${ROOT_LIBRARIES} ORStudyLib ${extra_libs}
    )
endfunction( SusyNtExec )

# build everything in util/ as an executable for ORStudy
file(GLOB files "util/*.cxx")
foreach(file ${files})
    SusyNtExec(${file})
endforeach()

atlas_install_data( data/* )
