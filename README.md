# Resilient Distributed Diffusion for Time-Dependent Byzantine Attacks
Simulation code of the paper ["Resilient Distributed Diffusion for Multi-task Estimation"](https://arxiv.org/abs/2003.11911) 
and ["Resilient Distributed Diffusion in Networks over Adversaries"](https://arxiv.org/pdf/2003.10563.pdf).

## Four scenarios are included as a combination of the following:
```
* Diffusion algorithm (DLMS) with stationary/non-stationary targets,
* Strong/weak attack as defined in "Resilient Distributed Diffusion in Networks over Adversaries",
* Small network for weak attack.
```
And all the four scenarios contain the attack, DLMSAW and Resilient-DLMSAW (R-DLMSAW). One can change the selection of compromised nodes and F values of F-DLMSAW.

We recommend to start with the case of strong attack with stationary targets. And proceed to other scenarios if interested.

## Directories
```
* "funcs" contains dependent functions,
* "plot_funcs" contains files to generate the figures in the paper, 
* "figures" contains .eps figures in the paper.
```

## Cite the paper
```
@ARTICLE{8926533,  
  author={J. {Li} and W. {Abbas} and X. {Koutsoukos}},  
  journal={IEEE Transactions on Signal and Information Processing over Networks},   
  title={Resilient Distributed Diffusion in Networks With Adversaries},   
  year={2020},  
  volume={6},    
  pages={1-17},}   

@INPROCEEDINGS{8510965,
  author={J. {Li} and X. {Koutsoukos}},  
  booktitle={2018 14th International Conference on Distributed Computing in Sensor Systems (DCOSS)},   
  title={Resilient Distributed Diffusion for Multi-task Estimation},   
  year={2018},    
  pages={93-102},}  
```
  
 
