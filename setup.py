"""
## Development

### release new version

    git commit -am "version bump";git push origin master
    python setup.py --version
    git tag -a v$(python setup.py --version) -m "upgrade";git push --tags

"""

import sys
if (sys.version_info[0]) != (3):
    raise RuntimeError('Python 3 required ')

import setuptools

with open('README.md', 'r') as fh:
    long_description = fh.read()

# install jupyter kernel
import sys
from setuptools.command.install import install

class PostInstallCommand(install):
    """
    Post-installation for installation mode.

    To automate the following:

        repo_name=$(python -c "import sys;print(sys.prefix.split('/')[-2])"); echo $repo_name;
        uv add --dev ipykernel;uv run ipython kernel install --user --name=$repo_name;
    """
    def run(self):
        install.run(self)
        kernel_name = "intom"  # Name of the kernel to install
        try:
            # List current kernels in JSON format
            output = subprocess.check_output(
                ["jupyter", "kernelspec", "list", "--json"],
                text=True
            )
            kernelspecs = json.loads(output).get("kernelspecs", {})
            if kernel_name in kernelspecs:
                print(f"Warning: Kernel '{kernel_name}' already exists. "
                      f"Kernel not installed. Use 'jupyter kernelspec list' to inspect existing kernels.")
            else:
                # Install the kernel if it does not exist
                subprocess.check_call([
                    sys.executable, "-m", "ipykernel", "install",
                    "--user", "--name", kernel_name, "--display-name", kernel_name
                ])
                print(f"Kernel '{kernel_name}' installed successfully.")
        except Exception as e:
            print(f"Error checking or installing kernel: {e}")

# g: Dynamically find subpackages in the 'modules' directory
_dynamic_packages = [f"intom.{_p}" for _p in setuptools.find_packages(where="modules")]
_dynamic_package_dir = {f"intom.{_p}": f"modules/{_p}" for _p in setuptools.find_packages(where="modules")}

# g: Dynamically find packages starting with 'scr_' in the root directory
_scr_dirs = setuptools.find_packages(where=".", include=["scr_*"])
_dynamic_scr_packages = [f"intom.{_p}" for _p in _scr_dirs]
_dynamic_scr_package_dir = {f"intom.{_p}": _p for _p in _scr_dirs}

setuptools.setup(
    long_description=long_description,
    long_description_content_type="text/markdown",
    cmdclass={'install': PostInstallCommand},

    # sync changes by uv pip install -e .
    packages=[
        'intom',
    ] + _dynamic_packages + _dynamic_scr_packages,     
    package_dir={
        'intom': 'modules',
        **_dynamic_package_dir,
        **_dynamic_scr_package_dir, # g: Unpacks the dynamically generated scr_ dictionary
    },
)