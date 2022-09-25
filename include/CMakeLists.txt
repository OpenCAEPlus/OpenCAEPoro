target_sources(
  OpenCAEPoro
  PUBLIC AllWells.hpp
         Doxygen.hpp
         FluidSolver.hpp
         LinearSystem.hpp
         Mixture.hpp
         OCPFluidMethod.hpp
         OCPStructure.hpp
         ParamOutput.hpp
         ParamWell.hpp
         UtilError.hpp
         Well.hpp
         BulkConn.hpp
         CornerGrid.hpp
         FaspSolver.hpp
         Grid.hpp
         MixtureBO.hpp
         OCPConst.hpp
         OCP.hpp
         OCPTable.hpp
         ParamRead.hpp
         Reservoir.hpp
         UtilInput.hpp
         WellPerf.hpp
         Bulk.hpp
         DenseMat.hpp
         FlowUnit.hpp
         LinearSolver.hpp
         MixtureComp.hpp
         OCPControl.hpp
         OCPOutput.hpp
         ParamControl.hpp
         ParamReservoir.hpp
         Solver.hpp
         UtilTiming.hpp)

target_include_directories(OpenCAEPoro PUBLIC ${CMAKE_CURRENT_LIST_DIR})
