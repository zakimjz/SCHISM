# SCHISM software : Authors - Karlton Sequeira, Mohammed Zaki (2004)

# INSTALLATION
there are three files: README, extract.sh and schism.tar.gz
To install, all files must be in the same directory say "/tmp/".
Run extract.sh while in this directory
It will unzip and untar schism.tar.gz and make the executables
It puts dataGen & convertData in "/tmp/schism/obj" & eclat in "/tmp/schism/src"

# FUNCTIONALITY
dataGen creates horizontal format ASCII high-dimensional datasets having
embedded subspaces

convertData converts horizontal format ASCII high-dimensional datasets output by
dataGen to a suitable format i.e. vertical/binary/IBM/WEKA/LDR/MAFIA having
embedded subspaces (Only IBM formats have been recently tested)

schism finds the embedded subspaces

The functionality of these algorithms are better described in the SCHISM paper:

K. Sequeira, M. J. Zaki. SCHISM: A New Approach for Interesting Subspace Mining.
In the Proceedings of the Fourth IEEE Conference On Data Mining. 2004.

# RUNNING THEM
Command line arguments for convertData and dataGen can be found by running them
without any arguments.

typically run them as

        dataGen -o /tmp/schism/data/swhatever.ha -d 60

This creates a horizontal ASCII dataset with default parameters and 60 dims.

        convertData -i /tmp/schism/data/swhatever -d 60

This creates a IBM format dataset with default parameters from the horizontal
ASCII file earlier created.

        schism -i /tmp/schism/data/swhatever.ibm

This mines the IBM format file for embedded subspaces using default parameters

It prints "<Execution time(minutes)> <# interesting subspaces found> <entropy of subspaces> <coverage>". To print the actual subspaces, add parameter " -o 1".

In order to run these programs using as few parameters involves suffixing the
horizontal ASCII file by ".ha" and placing all new files in the /tmp/schism/data
directory

Example use may be seen in the shell scripts dataSCHISM.sh and testSCHISM.sh. 
Paths must be appropriately changed to run the scripts.
