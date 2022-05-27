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
fi
if [ ! -d "venv/" ]; then
	echo "Python virtual environment will be created..."
	python3 -m venv --copies venv/
	source venv/bin/activate
	pip install -q -r requirements.txt
	echo "OK"
fi

echo "done!"
echo "don't forget to copy files into data/ and essential genes/"
