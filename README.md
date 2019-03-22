# Convenient Modeller
Convenient Modeller is a software which simplifies construction of kinetic models for cell metabolism, with optional integration of proteomics data. These are the source code.

# UPDATE:
Software is now available as an easy to use singularity image.

Singularity is a containerization software available for free here: https://singularity.lbl.gov/install-linux (instruction for linux)

After installation, users can download the image with the following line:

singularity pull shub://chuanfuyap/Convenient-Modeller

To run the image:

singularity run shub://chuanfuyap/Convenient-Modeller

# To install and run the software locally
This software requires the following libraries to run:-

Java Libraries:

Apache Common Maths version 3.5   http://archive.apache.org/dist/commons/math/

Java Native Access (JNA)          https://github.com/java-native-access/jna

Java Systems Biology Markup Language (JSBML) version 1    http://sbml.org/Software/JSBML



C Libraries: 
NOTE: it's possible the source file is corrupted and unable to be compiled, I would suggest deleting the functions in the source files that are giving errors during compilation step and compile again.
SOSLIB version 1.7 https://github.com/raim/SBML_odeSolver 

libSBML version 5.X.X https://sourceforge.net/projects/sbml/files/libsbml/5.13.0/ 

SUNDIALS version 2.3 http://computation.llnl.gov/projects/sundials/sundials-software

Relevant project:

NOTE: this C program must be compiled with the right library paths for GRaPe to be able to solve models in parameter estimation process. 
soslibJNIC https://github.com/chuanfuyap/soslibJNIC
