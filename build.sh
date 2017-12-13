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

echo "Compilation terminee."
