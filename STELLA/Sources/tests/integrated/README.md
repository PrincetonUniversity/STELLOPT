Integrated tests
================

The integrated tests here are written using the Python [pytest][pytest] package,
and make use of [xarray][xarray] for reading in the data. To install these
packages in a standalone environment, you can run:

    make create-test-virtualenv
    source tests/integrated/venv/bin/activate

This will create a Python [virtual environment][venv] with the packages needed
for running the tests. You can then run all the tests:

    make check-integrated-tests

This assumes that the executable is located two directories above this one. You
can set the environment variable `STELLA_EXE_PATH` to the location of the
`stella` executable if it has been built somewhere else.

Writing new tests
-----------------

The testing framework is setup to automatically find files called `test_*.py`
and run any functions it finds in them called `test_*`. Writing a new integrated
test for `stella` is as "simple" as writing a new function that starts with
`test_`.

Let's look at a basic test, `test_simple_run` in
[`simple/test_simple.py`](simple/test_simple.py) to see how to write a test that
runs `stella` and checks the output netCDF file:

```python
def test_simple_run(tmp_path):
    """Run a short test case to check that the output is generated with
    the expected fields, and the expected number of timesteps

    """
```

The function takes one argument called `tmp_path`: this is a path to a temporary
directory done magically by the test framework `pytest`. This is used so we can run
the test in an isolated environment, and not worry about other tests potentially
interfering with output files, for example.

We can help out future developers by providing a brief description of the aim of
the test, how it works, and so on, in the docstring.

The first three lines:

```python
    input_filename = "simple.in"
    copy_input_file(input_filename, tmp_path)
    os.chdir(tmp_path)
```

do some setting up by copying the input file we're using to the temporary
directory `pytest` has created for us, and then changing the working directory
to there. `copy_input_file` is a helper function defined in the same file.

Next we run `stella` with our input file using another helper function
`run_stella`:

```python
    run_stella(get_stella_path(), input_filename)
```

`run_stella` also checks that `stella` completed successfully. `get_stella_path`
is yet another helper function that returns the absolute path to the `stella`
executable.

We can now read the output file using `xarray`:

```python
    with xr.open_dataset(tmp_path / input_filename.replace("in", "out.nc")) as df:
```

`xarray` gives us a nice dictionary-style interface to the netCDF file, and has
some nice features like lazy loading -- only loading the datasets we actually
use, when we use them. Using the `with` statement like this just ensures that
the file gets closed properly should something unexpected happen in the rest of
the test. This can be important as netCDF files can be fussy about being closed.

The remainder of this function does the actual testing:

```python
        expected_keys = [
            "phi2",
            "phi_vs_t",
        ]

        for field in expected_keys:
            assert field in df

        assert len(df["t"]) == 11
        assert np.allclose(df["t"][-1], 10.0)
```

We use `assert` statements to check properties of the output file. In this test,
we check that the output file contains some expected diagnostics, that the
length of the time dimension is correct, and that the last timestep matches what
we expect.

[pytest]: https://pytest.org
[xarray]: http://xarray.pydata.org
[venv]: https://docs.python.org/3/library/venv.html
