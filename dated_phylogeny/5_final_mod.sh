# You have to change the following utility:

for file in ./*.xml;
do;
sed 's/beast.util.TreeParser/beast.base.evolution.tree.TreeParser/g' $file > $file.def;
mv $file.def $file;
done
