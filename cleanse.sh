#/bin/sh

mv data/initial_c.txt .;
mv data/bas.txt .;
mv data/atom.txt .;
mv data/env.txt .;
rm -rf data;
rm log;
mkdir data;
mv initial_c.txt data/initial_c.txt;
mv bas.txt data/bas.txt;
mv atom.txt data/atom.txt;
mv env.txt data/env.txt;
cp data/initial_c.txt data/para_c.txt;
echo 'Cleanse complete.You may restart your calculation.'
