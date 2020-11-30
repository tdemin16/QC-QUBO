g++ -std=c++2a -o3 -m64 -I/usr/include/python3.8 -fopenmp -Wall -Wextra main.cpp lib.cpp -o main

if [ $? -eq 0 ]; then
    echo ---Builded!---
    ./main
    echo ---Done---
else
    echo Build Failed!
fi 
