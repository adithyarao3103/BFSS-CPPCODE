include_directories(SYSTEM "${PROJECT_SOURCE_DIR}/Eigen/Eigen")

ADD_EXECUTABLE(testbfssconfigwriter testbfssconfigwriter.cpp)

ADD_EXECUTABLE(testbosonicaction testbosonicaction.cpp)

ADD_EXECUTABLE(bfsssimulation bfsssimulation.cpp)

ADD_EXECUTABLE(testdiracoperatorbfss testdiracoperatorbfss.cpp)

ADD_EXECUTABLE(testrationalapproximation testrationalapproximation.cpp)

ADD_EXECUTABLE(benchmark benchmark.cpp)