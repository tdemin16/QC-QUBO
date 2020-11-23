g++ -std=c++2a -o3 -m64 -I/usr/include/python3.8 -fopenmp -Wall -Wextra main.cpp lib.cpp -lpython3.8 -o main

if [ $? -eq 0 ]; then
    echo ---Builded!---
    ./main
    echo -----DONE DOING STUPID STUFF------
else
    echo Build Failed!
fi 
