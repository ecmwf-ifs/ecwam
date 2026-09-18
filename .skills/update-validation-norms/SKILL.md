---
name: Update validation norms
description: Instructions on how to update the reference validation norms for a specified set of tests and build precisions.
---

- Each ecWAM test yaml config in ecwam/tests/ contains a set of reference norms for the swh (significant wave height).
- From time-to-time, we may wish to update these references.
- The procedure for doing so is as follows:
    - Using the `package/bundle/arch/ecmwf/hpc2020/default` arch, do a clean `--build-type=debug` build.
    - Run the specified test(s). Start off with the unit tests configured to run via ctest in `tests/CMakeLists.txt`. Here, ignore the multi OMP thread and multi MPI rank tests for reference generation.
      Especially in debug builds, changing the launch config should not affect the output norms. For the remaining larger tests, use the contents of `cmake/ecwam_add_test.cmake` to fashion your own run
      scripts.
    - For each test, `<run-dir>/logs/model/statistics.log` contains the computed norms. Update the reference swh norms using computed values at the appropriate time.
    - Now run the tests again, and for the tests controlled via ctest, **do** run the multi OMP/MPI tests, and all computed norms should **bit-identically** match the reference.
    - Do a clean optimised build (simply remove the `--build-type` argument), and run the tests again (including multi OMP/MPI ctest tests), ensuring they all pass. This run is not
      expected to produce a bit-identical output.
    - Unless, otherwise specified, repeat this process for both DP and SP builds.
    - Preserve all the build dirs (at most there will be 4 accounting for the optimisation and precision combinations) so that the user
      can verify the new norms.
