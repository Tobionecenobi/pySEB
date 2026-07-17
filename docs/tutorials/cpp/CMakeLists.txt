set(PYSEB_DOC_TUTORIAL_NAMES
    diblock
    micelle
    nested-structures
    parameter-sweep
    pdb-polymer
    symbolic-expressions
    scattering-quantities
    figure13-star-chain
)

set(PYSEB_DOC_TUTORIAL_TARGETS)

foreach(tutorial_name IN LISTS PYSEB_DOC_TUTORIAL_NAMES)
    set(tutorial_target "pyseb-doc-${tutorial_name}")
    add_executable(
        ${tutorial_target}
        "${CMAKE_CURRENT_SOURCE_DIR}/${tutorial_name}.cpp"
    )
    target_compile_features(${tutorial_target} PRIVATE cxx_std_14)
    target_link_libraries(${tutorial_target} PRIVATE SEB::seb SEB::symbolic)
    set_target_properties(${tutorial_target} PROPERTIES
        OUTPUT_NAME "tutorial-${tutorial_name}"
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/docs-tutorials"
    )
    list(APPEND PYSEB_DOC_TUTORIAL_TARGETS ${tutorial_target})
endforeach()

add_custom_target(doc-tutorials DEPENDS ${PYSEB_DOC_TUTORIAL_TARGETS})
