# Resilient Distributed Diffusion for Time-Dependent Byzantine Attacks
## Simulation code of the paper ["Resilient Distributed Diffusion for Multi-task Estimation"](https://ieeexplore.ieee.org/document/8510965) and "Resilient Distributed Diffusion in Networks over Adversaries"

Four scenarios are included as a combination of the following:
---
* Diffusion algorithm (DLMS) with stationary/non-stationary targets,
* Strong/weak attack as defined in "Resilient Distributed Diffusion in Networks over Adversaries",
* Small network for weak attack.

And all the four scenarios contain the attack, DLMSAW and Resilient-DLMSAW (R-DLMSAW). One can change the selection of compromised nodes and F values of F-DLMSAW.

We recommend to start with the case of strong attack with stationary targets. And proceed to other scenarios if interested.

Directories:
---
* "funcs" contains dependent functions,
* "plot_funcs" contains files to generate the figures in the paper, 
* "figures" contains .eps figures in the paper.
