# Usage

## 1. Setup the venv
The MD_Engine uses a Python script to parse input, so you have to create a venv with the necessary dependencies. This only has to be done once per install.
```
cd MD_Engine
python3 -m venv venv
source venv/bin/activate
pip install -r src/py/requirements.txt
deactivate
```

## 2. Example run
```
make
./md_run input/2FBD.pdb
./md_run input/random_particles-1024.pdb
```# cpu_implementation
