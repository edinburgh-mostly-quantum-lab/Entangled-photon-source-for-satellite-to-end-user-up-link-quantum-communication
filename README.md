The scripts can be run either locally on your machine or using GitHub codespaces

# GitHub Codespaces (Requires GitHub account)
On this repo's GitHub page, click the green `Code` button, then the `codespaces` tab. Now click the green `Create codespace on master` button. This will create an environment hosted on GitHub with everything needed to run the scripts. The first launch may take some time as it installs texlive in the environment, which is a rather large dependency.

# Local
To run the scripts, first setup a python environment\
`python3 -m venv .venv`\
`source .venv/bin/activate`\
`pip install -r requirements.txt`

The scripts default to using tex, if tex is not installed on your system, then edit the line `plt.rc('text', usetex=True)` to `plt.rc('text', usetex=False)` in each script

Each script can then be run in with `python *.py` to produce a corresponding .png file\
The scripts can also be run in an interactive python session
