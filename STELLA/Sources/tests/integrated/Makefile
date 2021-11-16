# See README.md for more information on these integrated tests

# This directory
integrated_test_dir = $(GK_HEAD_DIR)/tests/integrated

# Create a virtualenv and install the test dependencies there
create-test-virtualenv:
	bash -c "cd $(integrated_test_dir); \
		python3 -m venv venv; \
		source venv/bin/activate; \
		pip install -r requirements.txt; \
		deactivate; \
	"

# Run the integrated tests
check-integrated: check-integrated-dependencies stella
	@echo "Running: integrated tests"
	@pytest

# Check that all the Python dependencies are importable
check-integrated-dependencies:
	@python3 $(integrated_test_dir)/check_integrated_requirements.py
