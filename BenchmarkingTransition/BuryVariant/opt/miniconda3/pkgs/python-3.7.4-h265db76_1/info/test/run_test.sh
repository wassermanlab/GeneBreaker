

set -ex



python -V
python3 -V
2to3 -h
pydoc -h
python3-config --help
python -c "import sysconfig; print(sysconfig.get_config_var('CC'))"
_CONDA_PYTHON_SYSCONFIGDATA_NAME=_sysconfigdata_x86_64_conda_cos6_linux_gnu python -c "import sysconfig; print(sysconfig.get_config_var('CC'))"
for f in ${PREFIX}/lib/python3.7/_sysconfig*.py; do echo "Checking $f:"; if [[ `rg @ $f` ]]; then exit 1; fi; done
exit 0
