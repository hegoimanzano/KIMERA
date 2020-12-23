# KIMERA

## Kinetic Monte Carlo code for Mineral Dissolution

KIMERA is a scientific tool for the study of mineral dissolution. It implements a reversible Kinetic Monte Carlo (KMC) method to study the time evolution of a dissolving system, obtaining the dissolution rate and information about the atomic scale dissolution mechanisms. KIMERA allows to define the dissolution process in multiple ways, using a wide diversity of event types to mimic the dissolution reactions, and define the mineral structure in great detail, including topographic defects, dislocations, and point defects. Therefore, KIMERA ensures to perform numerous studies with great versatility. In addition, it offers a good performance thanks to its parallelization and efficient algorithms within the KMC method. In this repository, we present the code features and show some examples of its capabilities. KIMERA is controllable via user commands, it is written in object-oriented C++, and it is distributed as open-source software.

KIMERA was written by Pablo Martin, and the v1.0 release was published in 2020.

To start dissolving your minerals with KIMERA, please check the github wiki: https://github.com/hegoimanzano/KIMERA/wiki. There, you will find instructions on how to compile the code and write your inpufile. There are also descriptions on the native KIMERA format used for restarting purposes, examples, and tests to guarantee the reproducibility.

Whenever you use KIMERA, please cite "Martin, P., Gaitero, J. J., Dolado, J. S., & Manzano, H. (2020). KIMERA: A Kinetic Montecarlo Code for Mineral Dissolution. Minerals, 10(9), 825" https://www.mdpi.com/2075-163X/10/9/825/htm
