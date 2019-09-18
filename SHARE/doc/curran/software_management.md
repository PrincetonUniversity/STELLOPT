Curran's guide to software management
=====================================

These suggestions strive to maximize traceability and minimize startup costs.  They are informed by experience working in two medium-sized collaborative software development organizations: Simulating eXtreme Spacetimes and SpaceX.  Implementing them requires discipline and some new infrastructure.

Ticket-driven development
-------------------------

All well-defined development initiatives should be tracked in tickets, preferably with relatively fine granularity (and development that is not clear in purpose should not be merged to master).  Every commit message must reference a ticket (enforced by a pre-merge gate).

Advantages of linking commits to tickets include:

* Visibility into new and ongoing efforts (avoid duplication, indicate active PoCs)
* Capture metadata such as related efforts and interested parties (by cross-linking tickets and commits and registering as watchers)
* Natural home for discussion related to the code (including tradeoffs that are unlikely to be captured in comments or commit messages)
* Natural home for tracking defects related to a merge (helping developers see where they made mistakes and to reinforce expectations)
* All code can be traced to a purpose, even if the original author is unavailable

Popular ticket systems include Jira, Trac, and GitHub Issues.

Reliable build system
---------------------

The build system (make, etc.) must be trustworthy, and ideally it should be performant as well.  Code should not be rebuilt if nothing changed, and it should always be rebuilt if something did.  Bit-for-bit reproducible builds and sandboxing build systems like Bazel are the ideal; parallel builds should not cause problems.  Running with stale or partial builds should be prevented (this typically requires that the build system be involved in running/deployment).

Build procedures should be well-documented and kept up-to-date, as should dependencies and other assumptions about the environment.  Build rules for a set of standard development/deployment environments should be maintained; standard containers can be used to minimize variance (and reduce maintenance burden).

A new user should be able to follow a single page of instructions (on any supported development platform) and be guaranteed to end up with a demonstrably working build.  Maintaining reliable supported environment definitions also facilitates continuous integration testing (below).

Golden master and continuous integration
----------------------------------------

The entire contents of the master branch should be maintained in an always-working state.  Whenever something new is merged to master, it becomes the responsibility of the entire collaboration to maintain it as other code is refactored.  In order to know that all code in master is working, it must be accompanied by some form of regression tests (spanning unit tests to integration tests), and those tests must be run regularly, flagging any merges that introduce regressions.

Code that does not build/run cleanly on master should either be deleted or moved to a branch.

Performance regression tests are extremely useful if the testing infrastructure has reliable performance characteristics (e.g. not an overloaded VM).

Documentation hygene
--------------------

Published documentation should be kept up-to-date or deleted (being version-controlled should go without saying).  Documentation should always build/render correctly.

Code should be documented using the language's conventions, and these unit-level docs should be made available and searchable on the web (e.g. Doxygen).  The web docs platform should be pared down to what is useful (that is, disable fancy features that interfere with the performance or reliability of delivering searchable author-provided docs).
