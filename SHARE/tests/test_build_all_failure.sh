#!/bin/sh
# Independent observable exit-status oracle; never compiles the full suite.
set -eu
solver_repo=$(CDPATH= cd -- "$(dirname "$0")/../.." && pwd)
task_scratch=$(mktemp -d /var/tmp/stellopt-build-exit-test.XXXXXXXX)
trap 'rm -rf "$task_scratch"' EXIT HUP INT TERM
mkdir -p "$task_scratch/bin" "$task_scratch/LIBSTELL" "$task_scratch/VMEC2000"
cat > "$task_scratch/bin/make" <<'MAKE'
#!/bin/sh
case "$FAIL_STAGE:$PWD:$*" in
    root:*:release) case "$PWD" in */LIBSTELL|*/VMEC2000) ;; *) exit 17 ;; esac ;;
    metadata:*:test_make) exit 19 ;;
    library:*/LIBSTELL:*) exit 23 ;;
    vmec:*/VMEC2000:*) exit 29 ;;
esac
exit 0
MAKE
chmod +x "$task_scratch/bin/make"
cd "$task_scratch"
for pair in none:0 root:17 metadata:19 library:23 vmec:29
do
    FAIL_STAGE=${pair%:*}
    expected=${pair#*:}
    export FAIL_STAGE
    if PATH="$task_scratch/bin:$PATH" "$solver_repo/build_all" \
        -o release -j 2 LIBSTELL VMEC2000 > "$task_scratch/log" 2>&1
    then
        observed=0
    else
        observed=$?
    fi
    if [ "$observed" != "$expected" ]
    then
        cat "$task_scratch/log"
        echo "FAIL: $FAIL_STAGE expected $expected, observed $observed" >&2
        exit 1
    fi
done
# Preserve normal skipping of unrequested codes; a blanket set -e would break it.
FAIL_STAGE=none PATH="$task_scratch/bin:$PATH" "$solver_repo/build_all" \
    -o release VMEC2000 > "$task_scratch/log" 2>&1
echo 'PASS: build_all propagates four failure stages and preserves selection'
