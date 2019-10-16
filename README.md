# SuicideBurnAutopilot

"Suicide Burn Autopilot" is a Python script using the Kerbal Space Program mods KRPC and MechJeb to land a spacecraft by performing a suicide burn maneuver including horizontal components. As no approximations are done in the calculations, the suicide burn is performed (most often) very accurately.


## Requirements

- Kerbal Space Program (KSP)
- KSP mods:
  - MechJeb2
  - KRPC
- Python 2.7
- Python packages:
  - numpy
  - scipy
  - krpc

Porting the script to Python 3.x should be comparably easy as only `print` and `xrange` have to be replaced. However, I did not test that.


## Description

A rocket burning retrograde in a gravitational field in vacuum is described by the equation of motion (EOM)

```
m(t) d²/dt² r(t) = -F*v_rot/|v_rot| - mu*m(t)*r(t)/|r(t)|³
```

where m(t) is the mass of the rocket, r(t) is the position vector of the rocket in the fixed body frame of the planet, F is the thrust of the rocket, v_rot is the velocity vector in the rotating body frame of the planet, and mu is the gravitational parameter of the planet. The script evaluates the flight trajectory of the suborbital rocket by numerically integrating the EOM and therefore calculating r(t). For the first part of the trajectory the rocket is coasting (F=0). At a certain time t_burn the rocket ignites the engine pointing retrograde and starts to slow down. A minimization algorithm then optimizes t_burn, so when the velocity reaches v_rot=0 at the end of the burn the altitude also is zero.

After calculating t_burn the script uses the KRPC autopilot to point retrograde and start the burn at the right time. The last few meters before touchdown are handled by the MechJeb autopilot.


## Usage

- Start KSP and fly a spacecraft.
- Make sure the KRPC server is running.
- Set RPC_PORT and STREAM_PORT to the right values (lines 32 and 33).
- F5
- Start the python script (usually `python suicide_burn_autopilot.py`).
- A small UI should appear on the right side of the screen.
- Deorbit your spacecraft, i.e. the periapsis has to be below the surface. Too shallow landing angles could lead to larger errors resulting in RUD events... but most of the time it still works pretty good.
- Hit the button "Start autopilot"
- The calculations will take about 30 to 60 seconds. You can pause the game in the meantime. Resume the game when the calculations are finished (indicated by the console output or the UI).
- Pray! The autopilot will do the rest.

In line 32 you can find the parameter FINAL_ALTITUDE which defines an offset above the surface in meters. The script will try to stop the rocket at that point and the MechJeb autopilot will then cover the last few meters. Increase FINAL_ALTITUDE for large rockets.

You can also try to mess with the integrating function by changing MAX_STEP_IVP. Generally, smaller values should give more precise results but will also take longer. However, the default values give a good compromise of computation time and not crashing into a planet.


## ~~Bugs~~Features

- The KRPC autopilot should point the rocket retrograde right after resuming the game (lines 375-377), but sometimes this is skipped. Does anybody know why?
- Why is the UI gone after landing but the Python process is still running?

If you find any new features or can help me to resolve some, please contact me. :)
