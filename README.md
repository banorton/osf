# Optical Simulation Framework
An Optical Simulation Framework (OSF) written in MATLAB for simulating simple systems using Fourier optics. It enables the simulation of light propagation through optical systems using various propagation methods, including the **Rayleigh-Sommerfeld integral** and the **Angular Spectrum method**. Users can define and manage optical elements like lenses, diffusers, and apertures, and simulate their effects on light fields in both 1D and 2D.

## Features
- **Simulation of Light Propagation**: Supports **Rayleigh-Sommerfeld** and **Angular Spectrum** methods for 1D and 2D fields.
- **Modular Optical Elements**: Add lenses, diffusers, and apertures to the optical system.
- **Customizable Fields**: Define fields with user-specified amplitude, phase, resolution, and dimensionality.
- **Advanced Visualization**:
  - Display amplitude and phase of fields.
  - Calculate and visualize Wigner Distribution Function (WDF).
  - Fourier transform visualization with `fft` and `dispFFT`.
- **Padding Support**: Adds padding for accurate field propagation.

## System Design
### 1. **`Sim` Class**: Simulation Controller
- Manages optical elements, field propagation, and system evolution.
- **Constructor**:
  - Initialize with parameters like resolution, field size, and wavelength.
- **Key Methods**:
  - `addLens`, `addAperture`, `addDiffuser`: Add elements to the system.
  - `prop`, `propToIndex`, `propToElement`, `propToDist`: Propagate light through the system.
### 2. **`Field` Class**: Light Field Representation
- Represents the light field with amplitude and phase.
- **Key Methods**:
  - `applyPhaseShift`, `applyPhaseRect`: Apply phase modifications.
  - `fft`, `dispFFT`: Compute and visualize Fourier transforms.
  - `wdf`, `dispWDF`: Calculate and display the Wigner Distribution Function.
### 3. **`Element` Class**: Abstract Optical Component
- Parent class for all optical elements.
- **Subclasses**:
  - `Lens`: Applies lens-based phase shifts.
  - `Diffuser`: Adds roughness and randomness to the phase.
  - `Aperture`: Applies spatial constraints to the field.

## Example: Implementation of a simple 4f system.
```matlab
resolution = 1e-6; fieldLength = 1e-3;
sim = Sim(resolution, fieldLength);

% Add lens for 4f system.
distanceAfterPreviousElement = 10e-3; lensFocalLength = 10e-3;
sim.addLens(distanceAfterPreviousElement, lensFocalLength, 'name', 'Lens 1');
sim.addLens(20e-3, 10e-3, 'name', 'Lens 2');

% Create a blank field with the correct dimensions for the simulation, then
% add a rectangular phase shift with a length of .2mm and a phase shift of pi.
field = sim.newField().applyPhaseRect(.2e-3, pi);

% Propagate through all elements and then 10mm past the last element.
% The optional verbose argument makes the function display the distances and
% elements at each step.
distAfterLastElement = 10e-3;
sim.prop(field, distAfterLastElement, 'verbose', true);
```

