g++ -std=c++2a -o3 -m64 -fopenmp -Wall -Wextra -D SIMULATION ./src/main.cpp ./src/lib.cpp ./src/generators.cpp -o ./bin/main

if [ ! -d "./bin" ]; then \
	mkdir ./bin; fi

if [ $? -eq 0 ]; then
    echo ---Builded!---
    cd bin/
    ./main
    echo -----END------
else
    echo Build Failed!
fi 