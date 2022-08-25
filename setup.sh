#!/bin/bash
cd "${0%/*}"
mkdir -p data
mkdir -p "essential genes"
if [ ! -d "blast/" ]; then
	echo "blast 2.12 will be downloaded..."
	wget -q "https://anaconda.org/bioconda/blast/2.12.0/download/linux-64/blast-2.12.0-hf3cf87c_4.tar.bz2"
		if [ $? -eq 0 ] 
		then 
			echo "OK" 
		else 
			echo "Could not download."
			exit 1
		fi
	mkdir blast
	tar xjf blast-2.12.0-hf3cf87c_4.tar.bz2 -C blast
	rm blast-2.12.0-hf3cf87c_4.tar.bz2
	if [[ ! $(ldconfig -p | grep libpcre.so.1) ]]; then
		# create missing symlink
		ln -s  $(find /lib/ -name libpcre.so.3 | head -n 1) blast/lib/ncbi-blast+/libpcre.so.1
	fi

fi
if [ ! -d "venv38/" ]; then
	echo "Python3 virtual environment will be created..."
	python3 -m venv --copies venv38/
	source venv38/bin/activate
	pip install -q -r requirements38.txt
	echo "OK"
fi


echo "done!"
echo "don't forget to copy files into data/ and essential genes/"
