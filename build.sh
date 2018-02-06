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
cp -f exec/DEMARENA testrelative/DEMARENA

#Move executable to public path
cp exec/DEMARENA ~/bin/
cp exec/POSTPROCESS ~/bin/

echo "Compilation terminee."
