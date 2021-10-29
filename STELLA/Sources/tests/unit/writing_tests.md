---
title: Writing Tests
---

[TOC]

This document will take you through writing a new [pFUnit][pfunit]
test, detailing some of the pitfalls and things to be aware of.

## Basic structure of a test

pFUnit can automatically detect and run new tests added to the test suite. In order to do
this, it uses a custom preprocessor to convert the `.pf` files to `.F90` which are then
compiled by the Fortran compiler as usual. We do also need to tell pFUnit about new files
however, and these are listed in `tests/pfunit_tests/testSuites.inc`.

A very simple pFUnit test for stella looks something like this:

`my_tests.pf`:

```
module test_basic_mod
  use PFUNIT_MOD
  implicit none

contains
  @test
  subroutine test_addition
    @assertEqual(4, 2 + 2)
  end subroutine test_addition
end module test_basic_mod
```

`testSuites.inc`:

```
ADD_TEST_SUITE(test_basic_mod_suite)
```

Note that the module name has an additional `_suite` appended in `testSuites.inc`.

Let's break down the test module line by line:

```
module test_basic_mod
```

The exact name isn't important, but in the stella tests, we prefer the naming convention
`test_<stella module name>_mod`; for example: `test_leq_mod`, or `test_stella_time_mod`.

```
  use PFUNIT_MOD
```

After the `.pf` to `.F90` conversion step, the test module will use various functions and
types from the pFUnit module, so we `use` the entire module.

The all-caps is important here: we currently support using either pFUnit 3.x or 4.x, which
use different module names. We preprocess the resulting `.F90` to replace `PFUNIT_MOD`
with the correct module name.

```
  implicit none
  
contains
```

In this basic test module we just have the one test procedure, but we could declare
module-level parameters here. For example, you may wish to set a `parameter` defining a
tolerance for equalities.

```
  @test
  subroutine test_addition
```

The `@test` is a pFUnit directive declaring the immediately following subroutine as a test
to be run automatically. The subroutine name is not important, but we prefer the naming
convention of a `test_` prefix to distinguish the procedure from non-test procedures in
the same module.

```
    @assertEqual(4, 2 + 2)
```

This is the basic pFUnit testing directive that checks for equality between its two
arguments. If they don't compare equal, pFUnit will flag the test as having failed and
report it at the end of the run. Other asserts include `@assertTrue`, `@assertFalse`,
and `@assertSameShape`.

Importantly, pFUnit will stop running _this particular_ test if an `assertEqual` fails. If
the arguments are arrays, pFUnit will report only the first index where they differ.

The order of the arguments doesn't particularly matter, but when reporting errors pFUnit
assumes the **first** argument is the _expected_ value, and the **second** argument is the
_actual_ value.

Something to be aware of here is that the pFUnit preprocessor doesn't like this directive
split over multiple lines! Currently, we just have to live with long lines.

Another thing to note is that all pFUnit directives (the keywords and functions beginning
with `@`) are case insensitive: `@AssertEqual` is the same as `@assertEqual` and
`@assertequal`.

## Getting fancier

### Setup and Teardown

After writing a few tests, you may find that there's repeated code setting things up
first, and then cleaning up afterwards. pFUnit allows you to declare subroutines to be run
before and after each test in a module:

```
@before
subroutine setup
  ! Do something before every test
end subroutine setup

@after
subroutine teardown
  ! Do something after every test
end subroutine teardown
```

The subroutine names are conventional, but `setup` and `teardown` are preferred.

**Important!** If you mark one subroutine with either `@before` or `@after`, you **must**
also have the corresponding one, even if it's just an empty routine.

These `@before` and `@after` routines are run before/after each individual test, even if
that test fails. This can be very useful to, for example, call the `<module>_finish`
routine to cleanup a stella module, which helps make the tests more independent. Without the
`@after` routine, a failure in one test might result in a stella module being in a bad or
wrong state for the next test. To guard against other tests or modules not cleaning up
after themselves correctly, you may wish to call the `<module>_finish` routine in both the
setup and teardown routines. This helps ensure the tests are independent.

### Tolerances

By default, pFUnit compares values exactly. This is usually not what you want when dealing
with floating point numbers. pFUnit can check equality within a tolerance, although
unfortunately, as of v4.1.9, only absolute tolerances (`abs(actual-expected) <=
tolerance`). Tolerances can be specified with the `tolerance` keyword to `@assertEqual`:

```
@assertEqual(3.14, pi, tolerance=1e-2)
```

### Message

When a test contains multiple asserts, it can be tricky to tell exactly which assert
caused a test to fail. Unfortunately, while pFUnit prints the values it's comparing, it
doesn't print the _expressions_. Instead, it has a mechanism to give a bit more
information with the `message` keyword:

```
@assertEqual(3.14, pi, message="comparing pi")
```

### MPI tests

pFUnit has the ability to _parameterise_ a test over the number of MPI ranks. That is, it
can rerun a test with different numbers of MPI processes. This can be useful to check that
a function works with different ranks.

Additionally, many routines often end up needing MPI to be set up because they contain
some code like `if (proc0) then ...`. For these routines to work correctly when called
from a pFUnit test, we essentially need to parameterise the test over one process.

The form of an MPI test is a little different:

```
@test(npes=[1, 2])
subroutine test_some_mpi(this)
  class(MPITestMethod), intent(inout) :: this
  ...
```

Here, we'll end up calling the test two times, first with one rank, and second with
two. Next, our test subroutine now takes a single argument of `class(MPITestMethod)` --
note it is `intent(inout)`. Conventionally we name it `this`, but the exact name is not
important. This argument has some useful methods:

```
this%getMpiCommunicator()    ! Get the MPI_COMM for this test run
this%getNumProcesses()       ! Wrapper for MPI_Comm_size
this%getProcessRank()        ! Wrapper for MPI_Comm_rank
```

Parameterised tests, including MPI tests, **must** have setup/teardown routines. The setup
routine **must** call `mp::init_mp`, but the teardown **must not** call `mp::finish_mp`:

```
@before
subroutine setup(this)
  use mp, only : init_mp
  class (MpiTestMethod), intent(inout) :: this
  ! Make sure stella variables like proc0 are populated correctly:
  call init_mp(this%getMpiCommunicator())
end subroutine setup

@after
subroutine teardown(this)
  class (MpiTestMethod), intent(inout) :: this
  ! No call to finish_mp!
end subroutine teardown
```

Calling `finish_mp` will finalise MPI, which will cause problems for the next MPI test!

### Parameterised tests

It's possible to parameterise tests over other parameters than just MPI ranks, though it
is a bit more involved and involves several pieces.

Here's an example from one of the stella parameterised tests:

```
module test_antenna_data_mod
  use PFUNIT_MOD

  implicit none

  @TestParameter
  !> Type to hold parameter we want to scan in
  type, extends(MpiTestParameter) :: test_antenna_parameter
    integer :: nk_stir_in
  contains
    procedure :: toString
  end type test_antenna_parameter

  @TestCase(constructor = new_case)
  !> Test fixture
  type, extends(MpiTestCase) :: test_antenna_case
    integer :: nk_stir_in
  contains
    procedure :: tearDown
  end type test_antenna_case

contains

  !> Construct test_antenna_case from test_antenna_parameter
  function new_case(test_parameter) result (test_case)
    type(test_antenna_case) :: test_case
    type(test_antenna_parameter), intent(in) :: test_parameter
    test_case%nk_stir_in = test_parameter%nk_stir_in
  end function new_case

  !> Construct test_antenna_parameter
  function new_parameter(nk_stir_in)
    type(test_antenna_parameter) :: new_parameter
    integer, intent(in) :: nk_stir_in
    new_parameter%nk_stir_in = nk_stir_in
  end function new_parameter

  !> Serialise test_antenna_parameter as a string
  function toString(this) result(string)
    class(test_antenna_parameter), intent(in) :: this
    character(:), allocatable :: string
    allocate(character(len=10)::string)
    write(string, '("nk_stir=", i2)') this%nk_stir_in
  end function toString

  !> Tear down test module
  subroutine tearDown(this)
    use antenna_data, only : finish_antenna_data
    class(test_antenna_case), intent(inout) :: this
    call finish_antenna_data
  end subroutine tearDown

  !> Define the parameters scan
  function get_parameters() result(params)
    type(test_antenna_parameter), allocatable :: params(:)

    params = [ &
         & new_parameter(nk_stir_in=-1), &
         & new_parameter(nk_stir_in=3) &
         & ]

  end function get_parameters

  @test(npes=[1,2,4], testParameters={get_parameters()})
  !> Test for initialization with no antenna
  subroutine test_init_antenna_data(this)
    use antenna_data
    implicit none
    class (test_antenna_case), intent(inout) :: this

    associate( nk_stir_in => this%nk_stir_in )
      call init_antenna_data(nk_stir_in)

      if (nk_stir_in < 0) then
        ! Antenna off
        @AssertFalse(ant_on)
      else
        ! Antenna on
        @AssertTrue(ant_on)
        @AssertEqual(nk_stir_in, nk_stir)
      end if
      @AssertSameShape(shape(nk_stir_in), shape(a_ant))
      @AssertSameShape(shape(nk_stir_in), shape(b_ant))
    end associate
  end subroutine test_init_antenna_data

end module test_antenna_data_mod
```

The parameterisation machinery looks very complicated here because we're just
parameterising over a single integer, with the result that we end up defining two very
similar types, but they are for different purposes, mostly to do with the internals of
pFUnit and how Fortran works.

Let's start with the test itself:

```
  @test(npes=[1,2,4], testParameters={get_parameters()})
  !> Test for initialization with no antenna
  subroutine test_init_antenna_data(this)
    use antenna_data
    implicit none
    class (test_antenna_case), intent(inout) :: this
```

Our `@test` directive now has parameterisation over both MPI ranks and an extra set of
parameters. The `testParameters` argument calls the function `get_parameters` to get an
array of _parameters_ which are used to construct _test cases_. The test cases are passed
into the test routine. pFUnit constructs the Cartesian product of `npes` and
`testParameters`, that is, it calls the test routine with every combination of `npes` and
`testParameters`.

**Note**: the curly brackets `{...}` around the argument to `testParameters` are required.

The function `get_parameters` looks like:

```
  function get_parameters() result(params)
    type(test_antenna_parameter), allocatable :: params(:)

    params = [ &
         & new_parameter(nk_stir_in=-1), &
         & new_parameter(nk_stir_in=3) &
         & ]

  end function get_parameters
```

This returns an `allocatable` array of `test_antenna_parameter` values. We want to scan
over the `integer` `nk_stir_in` -- why can't we just return an `integer` array? This is
because pFUnit needs to store the return value and has to know what type we're
using. Because Fortran doesn't have generic data types, we're forced to use a derived type
that extends pFUnit's `AbstractTestParameter`. Our definition of `test_antenna_parameter`:

```
  @TestParameter
  type, extends(MpiTestParameter) :: test_antenna_parameter
    integer :: nk_stir_in
  contains
    procedure :: toString
  end type test_antenna_parameter
```

actually extends `MpiTestParameter`, which itself extends `AbstractTestParameter`. We use
the `@TestParameter` directive to tell pFUnit this is our test parameter type. Notice
that we also define a type-bound procedure or method called `toString`. This is used to
print error messages. For our example, it can be quite simple:

```
  function toString(this) result(string)
    class(test_antenna_parameter), intent(in) :: this
    character(:), allocatable :: string
    allocate(character(len=10)::string)
    write(string, '("nk_stir=", i2)') this%nk_stir_in
  end function toString
```

The `MpiTestParameter` base class takes care of writing the number of processes, so we
just need to write the rest of the type to an `allocatable` `character`.

We also need to write a function to construct a `test_antenna_parameter` from an
`integer`:

```
  function new_parameter(nk_stir_in)
    type(test_antenna_parameter) :: new_parameter
    integer, intent(in) :: nk_stir_in
    new_parameter%nk_stir_in = nk_stir_in
  end function new_parameter
```

This isn't strictly necessary if you can make your array of test parameters some other
way, but this is usually the easiest and least error prone method.

Now that we have our test parameters, pFUnit uses these to construct a _test case_ which
is finally passed to the test routine. Our test case type looks like this:

```
  @TestCase(constructor = new_case)
  type, extends(MpiTestCase) :: test_antenna_case
    integer :: nk_stir_in
  contains
    procedure :: tearDown
  end type test_antenna_case
```

The `constructor` argument to the `@TestCase` directive tells pFUnit which function to use
to make a new `test_antenna_case`. For us, this is very similar to the constructor for
`test_antenna_parameter`:

```
  function new_case(test_parameter) result (test_case)
    type(test_antenna_case) :: test_case
    type(test_antenna_parameter), intent(in) :: test_parameter
    test_case%nk_stir_in = test_parameter%nk_stir_in
  end function new_case
```

We just copy the `nk_stir_in` member from the _parameter_ to the _case_.

Also notice that the `test_antenna_case` type has a `tearDown` method. Unfortunately for
parameterised tests, we can't just use the `@before` and `@after` directives. Instead, we
need to overload the `setUp` and `tearDown` methods on the test case type. We don't have a
`setUp` here, but we do have a `tearDown`, which just cleans up the module we're testing:

```
  subroutine tearDown(this)
    use antenna_data, only : finish_antenna_data
    class(test_antenna_case), intent(inout) :: this
    call finish_antenna_data
  end subroutine tearDown
```

Finally, now that we've constructed the parameters and the cases, we can actually use the
case in the test routine itself:

```
  subroutine test_init_antenna_data(this)
    use antenna_data
    implicit none
    class (test_antenna_case), intent(inout) :: this

    associate( nk_stir_in => this%nk_stir_in )
      call init_antenna_data(nk_stir_in)

      if (nk_stir_in < 0) then
        ! Antenna off
        @AssertFalse(ant_on)
      else
        ! Antenna on
        @AssertTrue(ant_on)
        @AssertEqual(nk_stir_in, nk_stir)
      end if
      @AssertSameShape(shape(nk_stir_in), shape(a_ant))
      @AssertSameShape(shape(nk_stir_in), shape(b_ant))
    end associate
  end subroutine test_init_antenna_data
```

We can access the `nk_stir_in` member of `test_antenna_case` through `this%nk_stir_in`. We
actually use the Fortran `associate` block here to make a local name, `nk_stir_in`, to
avoid writing `this%nk_stir_in` each time. The rest of the test looks like usual.

[pfunit]: https://github.com/Goddard-Fortran-Ecosystem/pFUnit
