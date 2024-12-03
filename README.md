To run the scripts, first setup a python environment\
`python3 -m venv .venv`\
`source .venv/bin/activate`\
`pip install -r requirements.txt`

The scripts default to using tex, if tex is not installed on your system, then edit the line `plt.rc('text', usetex=True)` to `plt.rc('text', usetex=False)` in each script

Each script can then be run in with `python *.py` to produce a corresponding .png file\
The scripts can also be run in an interactive python session
