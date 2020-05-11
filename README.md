# Convenient Modeller
Convenient Modeller is a software which simplifies construction of kinetic models for cell metabolism (with inclusion of regulatory effects), with optional integration of proteomics data. 

The tool allows for use of metabolite concentration, and/or flux for fitting a model, main objective of this tool is to allow for rapid prototyping of the a metabolic phenomenon of interest with a kinetic/dynamic model. For example, if one is interested in glycolysis of a species such as baker's yeast, metabolite/flux values can be obtained (from literature of experiment) to fit a model. Following the fit, a **[SBML]** file  is generated which can be exported to any SBML simulator such as **[COPASI]**, **[SBMLsimulator]** or **[JWS online]** to be manipulated (e.g. removal of enzymes, removal of activators/inhibitors or increasing enzyme concentration) and simulated to see what its outcome would be before proceeding with next step of metabolic study.

[SBML]: http://sbml.org/Main_Page
[COPASI]: http://copasi.o:
[SBMLsimulator]: https://github.com/draeger-lab/SBMLsimulator
[JWS online]: https://jjj.mib.ac.uk


# How to use Convenient Modeller
NOTE: the GUI has bugs that has yet to be patched, but the optimisation engine that runs it is bug free. GUI bugs are mostly errors messages but should not affect model building. If GUI bugs are preventing model construction, the TSV format (see below) is recommended for model building and running the 

## Option 1 (Recommended): Singularity
Software is now available as an easy to use singularity image.

Singularity is a containerization software available for free here: https://singularity.lbl.gov/install-linux (instruction for linux)

After installation, users can download the image with the following line:

singularity pull shub://chuanfuyap/Convenient-Modeller

To run the image:

singularity run shub://chuanfuyap/Convenient-Modeller

### Mac users:
In order for Mac users to use the GUI version of the software via singularity, they need to enable X11 Forwarding in vagrant box, as well as install xauth server within the vagrant box

##  Option 2: To install and run the software locally
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

NOTE: this C program must be compiled with the right library paths for Convenient-Modeller to be able to solve models in parameter estimation process. 
soslibJNIC https://github.com/chuanfuyap/soslibJNIC

# HOW TO USE:

Full tutorial available **[here]**
[here]: https://github.com/chuanfuyap/Convenient-Modeller/blob/master/Convenient%20Modeller%20Tutorial.pdf

## GUI
Crash Course

<img src="https://github.com/chuanfuyap/Convenient-Modeller/blob/master/images/crash_course.png"
  alt="Size Limit comment in pull request about bundle size changes"
  width="960" height="580">

## Alternative model building with TSV
An alternative method is to use a spreadsheet to build the model before transfering to a text file, which can then be exported the convenient modeller GUI for further manipulation or model fitting.
