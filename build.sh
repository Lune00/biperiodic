#Recompile tout le projet

#clen lib
rm lib/*.o

#Compile les sources:
echo "source"
cd src
make clean
make

cd ../
cd exec
echo "app"
make clean
make

cd ..

#Copie executable dans un dossier d'essai
mkdir -p test
#cp exec/BPDEM2D test/BPDM2D
cp -f exec/DEMARENA test/DEMARENA
cp -f exec/DEMARENA testcollisions/DEMARENA
cp -f exec/DEMARENA testshear3x3/DEMARENA
cp -f exec/DEMARENA testshear3x3_2/DEMARENA
cp -f exec/DEMARENA testcompression/DEMARENA
cp -f exec/DEMARENA testreal/DEMARENA
cp -f exec/DEMARENA testfriction/DEMARENA
cp -f exec/DEMARENA test_shearIsoCompressed/DEMARENA
cp -f exec/DEMARENA testcompression2/DEMARENA
cp -f exec/DEMARENA testcompression3/DEMARENA
cp -f exec/DEMARENA testcompression4/DEMARENA
cp -f exec/DEMARENA testcompression4_frictionless/DEMARENA
cp -f exec/DEMARENA testForce/DEMARENA
cp -f exec/DEMARENA testForce2/DEMARENA
cp -f exec/DEMARENA testForce2/drive/DEMARENA

echo "Compilation terminee."
