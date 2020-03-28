mkdir schism
cd schism
mkdir src obj data tmp
cd tmp
mv ../../schism.tar.gz .
tar -zPxvf schism.tar.gz
g++ -O4 -o ../obj/dataGen ../src/dataGen.cc
g++ -O4 -o ../obj/convertData ../src/convertData.cc
cd ../src
make
cd ../tmp
