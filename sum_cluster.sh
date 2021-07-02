
for i in seq -f '%04g' 0 99;
do
	cd $i
  cp * ../.
  cd ../
done
