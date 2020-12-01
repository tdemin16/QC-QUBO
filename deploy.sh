BUILD="g++ -std=c++2a -o3 -m64 -I/usr/include/python3.8 -fopenmp -Wall -Wextra ./src/main.cpp ./src/lib.cpp ./src/generators.cpp -o ./bin/main"

if [ ! -d "./bin" ]; then \
	mkdir ./bin; 
fi

$BUILD

if [ $? -eq 0 ]; then
    echo ---Builded!---
    cd bin/
    ./main
    echo ---Done---
else
    echo Build Failed!
fi 
