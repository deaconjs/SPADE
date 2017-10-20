### SPADE
#### Structural Proteomics Application Development Environment 

[https://sites.google.com/view/spade](sites.google.com/view/spade)


#### Linux:

[Here](https://danieleriksson.net/2017/02/08/how-to-install-latest-python-on-centos/) are instructions for setting up a good Python build on CentOS/RH. 

##### Install SPADE

SPADE requires build-essential, python-tk, python-vtk, python-pmw, python-numpy, and python-scipy. I installed all those on linux mint with apt-get. Then install biopython with "sudo pip install biopython". 

##### Launch SPADE 

$ python SPADE.py

##### Launch MolecularViewer:

$ python MolecularViewer.py

##### Install AlignmentEditor

Close SPADE then cd into the Applications folder. Clone the AlignmentEditor github repo into a folder named AlignmentEditor. Then start SPADE.


#### Windows:

Windows installer [here](https://sourceforge.net/projects/spade/). This version has limited functionality.
