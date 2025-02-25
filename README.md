# Package of Numerical Weather Prediction models: Bolam, Moloch, Globo.

# Version 24.1, October 2024: Change of format of all internal input-output
files (unmormatted binary): mhf, shf, model_param_constant.bin,
time consuming better for 10 %

Bolam and Globo models are models based on hydrostatic approach for the dynamics, while Moloch is based on a non-hydrostatic approach.

Bolam and Moloch are limited area models, and Globo covers the global domain.

All models use the Eulerian numerical approach implemented on regular lon-lat grid of rotated geographical coordinate system with Arakawa-C staggering.

Bolam and Globo use terrain-following sigma atmospheric vertical coordinate, Moloch uses a hybrid terrain-following coordinate «zita», relaxing smoothly to horizontal surfaces away from the earth surface.

These models include state-of-the-art parametrizations of physical processes, modules for controlling input/output and run configuration.

Models are suitable for operational weather forecast and for meteorological research studies.

The implementation instruction can be found [here](isac-nwp-models_readme.txt)
