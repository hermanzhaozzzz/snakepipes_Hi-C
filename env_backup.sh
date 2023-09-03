mamba env export > conda_env.yml
pip freeze | grep -v conda > pip_env.txt
