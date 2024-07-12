MODULE GIT_VERSION_MOD
   IMPLICIT NONE
#if defined(GIT_VERSION_EXT)
   CHARACTER(64), PARAMETER :: git_repository = GIT_REPO_EXT
   CHARACTER(32), PARAMETER :: git_version = GIT_VERSION_EXT
   CHARACTER(40), PARAMETER :: git_hash = GIT_HASH_EXT
   CHARACTER(32), PARAMETER :: git_branch = GIT_BRANCH_EXT
   CHARACTER(19), PARAMETER :: built_on = BUILT_ON_EXT
#else
   CHARACTER(64), PARAMETER :: git_repository = "not from a git repo"
   CHARACTER(32), PARAMETER :: git_version = ""
   CHARACTER(40), PARAMETER :: git_hash = ""
   CHARACTER(32), PARAMETER :: git_branch = ""
   CHARACTER(19), PARAMETER :: built_on = ""
#endif
END MODULE GIT_VERSION_MOD
