#!/usr/bin/env bash

# define the default test to run (defined in the tests directory)
default_test=etopo1_oper_an_fc_O48_cy50r1.yml
test=${1:-$default_test}

# run preprocessor (usually only has to be done once)
build/bin/ecwam-run-preproc  --config=tests/$test

# run presets (usually only has to be done once)
build/bin/ecwam-run-preset   --config=tests/$test

# run model
build/bin/ecwam-run-model    --config=tests/$test

# sweep
exp=$(basename $test .yml)
rundir=../ecwam_runs/${exp}
for dir in logs restart output; do
	mkdir -p "$rundir/$dir"
	cp -a "$dir"/* "$rundir/$dir/"
done

echo "Model run complete. Output is in $rundir"