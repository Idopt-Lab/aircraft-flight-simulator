1. Add the file to the MATLAB path.
2. Run F16Lookup.m to load the data.
3. Run simulinkscript.m
- If Simulink isn’t open, manually open Flightsim2.slx.
- Ensure the initial conditions are loaded in the 6DOF block.
- Note: The simulation may take up to 10 minutes to complete.
3. Open a terminal in the project folder (ensure Python is installed on your system).
4. Create a virtual environment:
- On Windows: python -m venv venv
- On Linux: python3 -m venv venv
5. Activate the virtual environment:
- On Windows: .\venv\Scripts\Activate
- On Linux: source .venv/bin/activate
6. Run python main.py.
7. Open http://localhost:8000 in any browser.
8. Rerun the file(Simulink Script)
