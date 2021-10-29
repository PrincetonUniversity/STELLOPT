# This script just checks if we can import all the modules listed in
# requirements.txt

import importlib
import pathlib
import textwrap

requirements_file = pathlib.Path(__file__).parent / "requirements.txt"

with open(requirements_file, "r") as f:
    contents = f.readlines()

requirements = []
for line in contents:
    if line.startswith("#"):
        continue
    module = line.split("=")[0]
    requirements.append(module.replace(">", "").replace("~", ""))

failed_modules = []
for requirement in requirements:
    try:
        importlib.import_module(requirement)
    except ImportError:
        failed_modules.append(requirement)

if failed_modules:
    list_of_failed_modules = ", ".join(failed_modules)
    print(
        textwrap.dedent(
            f"""\
            Could not import: {list_of_failed_modules}
            Run:
                make create-test-virtualenv
                source tests/integrated/venv/bin/activate
            to install and activate test dependencies
            """
        )
    )
    exit(1)
