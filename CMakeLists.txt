cmake_minimum_required(VERSION 3.28)

project(CoupledPhysicsSolver VERSION 0.1.0 LANGUAGES CXX)

if(NOT APPLE)
    message(FATAL_ERROR "CoupledPhysicsSolver requires Apple Metal on macOS")
endif()

enable_language(OBJCXX)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include(CTest)
find_library(FOUNDATION_FRAMEWORK Foundation REQUIRED)
find_library(METAL_FRAMEWORK Metal REQUIRED)
find_library(ACCELERATE_FRAMEWORK Accelerate REQUIRED)

add_subdirectory(matter)

set(COUPLED_SHADER_DIR "${CMAKE_CURRENT_BINARY_DIR}/shaders")
set(COUPLED_METALLIB "${COUPLED_SHADER_DIR}/CoupledPhysics.metallib")
set(COUPLED_AIR_FILES)

function(coupled_metal_air name source)
    cmake_parse_arguments(ARG "PRECISE" "" "DEFINES" ${ARGN})
    set(output "${COUPLED_SHADER_DIR}/${name}.air")
    set(options -std=metal4.0 -O3)
    if(ARG_PRECISE)
        list(APPEND options -fno-fast-math)
    endif()
    foreach(definition IN LISTS ARG_DEFINES)
        list(APPEND options "-D${definition}")
    endforeach()
    add_custom_command(
        OUTPUT "${output}"
        COMMAND "${CMAKE_COMMAND}" -E make_directory "${COUPLED_SHADER_DIR}"
        COMMAND xcrun -sdk macosx metal ${options}
            -I "${CMAKE_CURRENT_SOURCE_DIR}/include"
            -c "${CMAKE_CURRENT_SOURCE_DIR}/${source}"
            -o "${output}"
        DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/${source}"
        VERBATIM
    )
    set(COUPLED_AIR_FILES ${COUPLED_AIR_FILES} "${output}" PARENT_SCOPE)
endfunction()

coupled_metal_air(ConstraintSolver src/metal/ConstraintSolver.metal)
coupled_metal_air(QualityContact src/metal/MetalQualityContactSolver.metal PRECISE)
coupled_metal_air(UnifiedQuality src/metal/MetalUnifiedQualitySolver.metal PRECISE)
coupled_metal_air(FreeBody src/metal/FreeBodyDynamics.metal)
coupled_metal_air(Multicopter src/metal/MulticopterDynamics.metal)
coupled_metal_air(Collision src/metal/Collision.metal PRECISE)
coupled_metal_air(Broadphase src/metal/DeterministicBroadphase.metal)
coupled_metal_air(ArticulatedOperator src/metal/ArticulatedOperator.metal)
coupled_metal_air(ParameterizedArticulatedOperator src/metal/ArticulatedOperator.metal
    DEFINES MR_ARTICULATED_OPERATOR_KERNEL_NAME=mr_parameterized_articulated_operator MR_ARTICULATED_OPERATOR_BODY_PARAMETERS=1)
coupled_metal_air(ArticulatedABA src/metal/ArticulatedABA.metal)
coupled_metal_air(ParameterizedArticulatedABA src/metal/ArticulatedABA.metal
    DEFINES MR_ABA_KERNEL_NAME=mr_parameterized_articulated_aba_step MR_ABA_BODY_PARAMETERS=1)
coupled_metal_air(ArticulatedABASmall src/metal/ArticulatedABA.metal
    DEFINES MR_ABA_KERNEL_NAME=mr_articulated_aba_step_small MR_ABA_COMPILED_MAX_BODIES=12 MR_ABA_COMPILED_MAX_DOFS=16 MR_ABA_COMPILED_MAX_Q=17)
coupled_metal_air(MultiArticulatedABA src/metal/ArticulatedABA.metal
    DEFINES MR_ABA_KERNEL_NAME=mr_multi_articulated_aba_step MR_ABA_MULTI_ARTICULATION=1)
coupled_metal_air(ParallelArticulatedABA src/metal/ArticulatedParallelABA.metal)
coupled_metal_air(ParallelParameterizedArticulatedABA src/metal/ArticulatedParallelABA.metal
    DEFINES MR_PARALLEL_ABA_KERNEL_NAME=mr_parallel_parameterized_articulated_aba_step MR_ABA_BODY_PARAMETERS=1)
coupled_metal_air(ParallelArticulatedABASmall src/metal/ArticulatedParallelABA.metal
    DEFINES MR_PARALLEL_ABA_KERNEL_NAME=mr_parallel_articulated_aba_step_small MR_ABA_COMPILED_MAX_BODIES=12 MR_ABA_COMPILED_MAX_DOFS=16 MR_ABA_COMPILED_MAX_Q=17)
coupled_metal_air(ParallelMultiArticulatedABA src/metal/ArticulatedParallelABA.metal
    DEFINES MR_PARALLEL_ABA_KERNEL_NAME=mr_parallel_multi_articulated_aba_step MR_ABA_MULTI_ARTICULATION=1)
coupled_metal_air(DiscreteElasticRod src/metal/DiscreteElasticRod.metal PRECISE)
coupled_metal_air(GeneralizedConstraints src/metal/MetalGeneralizedConstraints.metal PRECISE)
coupled_metal_air(MultiArticulatedContact src/metal/MetalMultiArticulatedContact.metal PRECISE)
coupled_metal_air(ArticulatedInverseMass src/metal/ArticulatedInverseMass.metal)
coupled_metal_air(MultiArticulatedInverseMass src/metal/ArticulatedInverseMass.metal
    DEFINES MR_INVERSE_MASS_KERNEL_NAME=mr_multi_articulated_inverse_mass MR_INVERSE_MASS_MULTI_ARTICULATION=1)
coupled_metal_air(ParameterizedMultiArticulatedInverseMass src/metal/ArticulatedInverseMass.metal
    DEFINES MR_INVERSE_MASS_KERNEL_NAME=mr_parameterized_multi_articulated_inverse_mass MR_PARALLEL_INVERSE_MASS_KERNEL_NAME=mr_parameterized_parallel_multi_articulated_inverse_mass MR_INVERSE_MASS_MULTI_ARTICULATION=1 MR_INVERSE_MASS_BODY_PARAMETERS=1)
coupled_metal_air(StreamingArticulatedInverseMass src/metal/ArticulatedInverseMass.metal
    DEFINES MR_INVERSE_MASS_KERNEL_NAME=mr_world_streaming_articulated_inverse_mass MR_PARALLEL_INVERSE_MASS_KERNEL_NAME=mr_world_parallel_streaming_articulated_inverse_mass MR_INVERSE_MASS_BODY_PARAMETERS=1 MR_INVERSE_MASS_STREAMING_RHS=1 MR_PARALLEL_INVERSE_MASS_STREAMING_RHS=1)
coupled_metal_air(World src/metal/MetalWorld.metal)
coupled_metal_air(WorldContact src/metal/MetalWorldContact.metal)
coupled_metal_air(WorldCompiler src/metal/WorldCompiler.metal)
coupled_metal_air(LocomotionTask src/metal/LocomotionTask.metal)
coupled_metal_air(PolicyInference src/metal/PolicyInference.metal)
coupled_metal_air(Metal4ExecutionProbe src/metal/Metal4ExecutionProbe.metal PRECISE)

add_custom_command(
    OUTPUT "${COUPLED_METALLIB}"
    COMMAND xcrun -sdk macosx metallib ${COUPLED_AIR_FILES} -o "${COUPLED_METALLIB}"
    DEPENDS ${COUPLED_AIR_FILES}
    VERBATIM
)
add_custom_target(coupled_shaders ALL DEPENDS "${COUPLED_METALLIB}")

add_library(coupled_physics STATIC
    src/evidence/TrialEvidence.mm
    src/metal/Metal4Execution.mm
    src/core/Model.cpp
    src/core/CpuDynamics.cpp
    src/core/EngineModel.cpp
    src/core/EngineModelComposer.cpp
    src/core/Kinematics.cpp
    src/core/ParallelABASchedule.cpp
    src/core/CompiledMetalMultiArticulatedProgram.cpp
    src/core/CompiledMetalMultiArticulatedContactProgram.cpp
    src/core/WorldCompiler.cpp
    src/core/WorldPack.cpp
    src/core/GeometryCooker.cpp
    src/core/SurgicalPSM.cpp
    src/core/Wave.cpp
    src/core/SurgicalAssets.cpp
    src/core/SurgicalWorld.cpp
    src/core/HeterogeneousWorld.cpp
    src/core/DiscreteElasticRod.cpp
    src/core/Collision.cpp
    src/core/FreeBodyDynamics.cpp
    src/core/Multicopter.cpp
    src/core/ArticulatedDynamics.cpp
    src/core/ArticulatedActuation.cpp
    src/core/ArticulatedJointLimits.cpp
    src/core/ArticulatedContact.cpp
    src/core/CoupledArticulatedRigidContact.cpp
    src/core/ArticulatedRigidCollision.cpp
    src/core/ArticulatedRigidWorld.cpp
    src/core/ArticulatedMixedConstraints.cpp
    src/core/ArticulatedCollision.cpp
    src/core/ArticulatedWorld.cpp
    src/core/MultiArticulatedWorld.cpp
    src/core/MultiArticulatedContact.cpp
    src/core/ConstraintIR.cpp
    src/core/ReferenceConicSolver.cpp
    src/core/QualityContactSolver.cpp
    src/core/ConstraintSolver.cpp
    src/core/RigidBodyWorld.cpp
    src/core/TaskProgram.cpp
    src/core/PolicyProgram.cpp
    src/core/InteractionPack.cpp
    src/metal/MetalArticulatedOperator.mm
    src/metal/MetalArticulatedABA.mm
    src/metal/MetalMultiArticulatedABA.mm
    src/metal/MetalDiscreteElasticRod.mm
    src/metal/MetalArticulatedInverseMass.mm
    src/metal/MetalMultiArticulatedInverseMass.mm
    src/metal/MetalMultiArticulatedConstraints.mm
    src/metal/MetalMultiArticulatedContact.mm
    src/metal/MetalQualityContactSolver.mm
    src/metal/MetalUnifiedQualitySolver.mm
    src/metal/MetalWorld.mm
    $<TARGET_OBJECTS:numi_matter_compiler_objects>
    $<TARGET_OBJECTS:numi_matter_runtime_objects>
)
add_dependencies(coupled_physics coupled_shaders numi_matter_shaders)
target_include_directories(coupled_physics PUBLIC include matter/include PRIVATE src)
target_compile_definitions(coupled_physics PRIVATE
    METALROBO_BUILDING_LIBRARY=1
    METALROBO_DEFAULT_METALLIB="${COUPLED_METALLIB}")
target_compile_options(coupled_physics PRIVATE
    "$<$<COMPILE_LANGUAGE:CXX>:-Wall;-Wextra;-Wpedantic>"
    "$<$<COMPILE_LANGUAGE:OBJCXX>:-fobjc-arc;-Wall;-Wextra;-Wpedantic>")
target_link_libraries(coupled_physics PUBLIC
    "${FOUNDATION_FRAMEWORK}" "${METAL_FRAMEWORK}" "${ACCELERATE_FRAMEWORK}")

function(coupled_probe target source)
    add_executable(${target} ${source})
    add_dependencies(${target} coupled_shaders numi_matter_shaders)
    target_link_libraries(${target} PRIVATE coupled_physics)
    target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic -Werror)
    set_target_properties(${target} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")
endfunction()

coupled_probe(coupled_matter_physics_probe apps/matter_physics_probe.mm)
target_compile_options(coupled_matter_physics_probe PRIVATE -fobjc-arc)
target_compile_definitions(coupled_matter_physics_probe PRIVATE
    NUMI_MATTER_MATERIAL="${CMAKE_CURRENT_SOURCE_DIR}/matter/materials/silicone.nmatter"
    NUMI_MATTER_STATEFUL_MATERIAL="${CMAKE_CURRENT_SOURCE_DIR}/matter/materials/damageable_silicone.nmatter"
    NUMI_MATTER_METALLIB="${CMAKE_CURRENT_BINARY_DIR}/matter/shaders/NumiMatter.metallib")

coupled_probe(coupled_surgical_tissue_probe apps/surgical_tissue_probe.mm)
target_compile_options(coupled_surgical_tissue_probe PRIVATE -fobjc-arc)
target_compile_definitions(coupled_surgical_tissue_probe PRIVATE
    NUMI_JEJUNUM_MATERIAL="${CMAKE_CURRENT_SOURCE_DIR}/matter/materials/porcine_jejunum_fung.nmatter"
    NUMI_MATTER_METALLIB="${CMAKE_CURRENT_BINARY_DIR}/matter/shaders/NumiMatter.metallib")

coupled_probe(coupled_suture_handoff_probe apps/dual_psm_suture_handoff_probe.cpp)
target_compile_definitions(coupled_suture_handoff_probe PRIVATE
    NUMI_JEJUNUM_MATERIAL="${CMAKE_CURRENT_SOURCE_DIR}/matter/materials/porcine_jejunum_fung.nmatter"
    NUMI_MATTER_METALLIB="${CMAKE_CURRENT_BINARY_DIR}/matter/shaders/NumiMatter.metallib")

coupled_probe(coupled_trial_evidence_probe apps/trial_evidence_probe.mm)
target_compile_options(coupled_trial_evidence_probe PRIVATE -fobjc-arc)

coupled_probe(coupled_trial_check apps/trial_check.cpp)

coupled_probe(
    coupled_perfused_active_tissue_contract_probe
    apps/perfused_active_tissue_contract_probe.cpp)
target_compile_definitions(
    coupled_perfused_active_tissue_contract_probe
    PRIVATE NUMI_JEJUNUM_MATERIAL="${CMAKE_CURRENT_SOURCE_DIR}/matter/materials/porcine_jejunum_fung.nmatter")

coupled_probe(coupled_metal4_execution_probe apps/metal4_execution_probe.mm)
target_compile_options(coupled_metal4_execution_probe PRIVATE -fobjc-arc)
target_compile_definitions(coupled_metal4_execution_probe PRIVATE
    COUPLED_METALLIB="${COUPLED_METALLIB}")

if(BUILD_TESTING)
    add_test(NAME coupled.matter.metal_world
        COMMAND $<TARGET_FILE:coupled_matter_physics_probe> --metal-world-coupling)
    add_test(NAME coupled.suture.tissue_coupling
        COMMAND $<TARGET_FILE:coupled_suture_handoff_probe> --tissue-coupling-only)
    add_test(NAME coupled.evidence.trial_contract
        COMMAND $<TARGET_FILE:coupled_trial_evidence_probe>)
    add_test(NAME coupled.tissue.perfused_active_contract
        COMMAND $<TARGET_FILE:coupled_perfused_active_tissue_contract_probe>)
    add_test(NAME coupled.metal4.execution_contract
        COMMAND $<TARGET_FILE:coupled_metal4_execution_probe>)
    set_tests_properties(
        coupled.matter.metal_world
        coupled.suture.tissue_coupling
        coupled.evidence.trial_contract
        coupled.tissue.perfused_active_contract
        coupled.metal4.execution_contract
        PROPERTIES LABELS "coupled;metal;physics" TIMEOUT 300)
endif()
