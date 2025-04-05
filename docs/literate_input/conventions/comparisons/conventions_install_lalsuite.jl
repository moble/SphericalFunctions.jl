# Construct the CondaPkg.toml file to use to make sure we get the right Python version and
# we get pip installed.
conda_pkg_toml = """
[deps]
python = "<3.13"
pip = ""
numpy = "==2.2.4"
"""
try
    open(joinpath(LOAD_PATH[2], "CondaPkg.toml"), "w") do io
        write(io, conda_pkg_toml)
    end
catch e
    println("Error copying CondaPkg.toml: $e")
end

# Now we'll set up the CondaPkg environment
import CondaPkg
import PythonCall

# This ugly hack is to ensure that lalsuite is installed without any dependencies; by
# default it comes with lots of things we don't need that break all the time, so I really
# don't want to bother fixing them.  The `--no-deps` flag is not supported by CondaPkg, so
# we have to use PythonCall to install it.
PythonCall.@pyexec `
from sys import executable as python;
import subprocess;
subprocess.call([python, "-m", "pip", "install", "-q", "--no-deps", "lalsuite==7.25.1"]);
`
