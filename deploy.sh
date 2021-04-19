BUILD="g++ -std=c++2a -o3 -m64 -g -fsanitize=address -fopenmp -Wall -Wextra ./src/main.cpp ./src/lib.cpp ./src/npp.cpp ./src/qap.cpp ./src/tsp.cpp -o ./bin/main"

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
