g++ -std=c++2a -o3 -fopenmp -Wall -Wextra -D SIMULATION -D DEBUG main.cpp lib.cpp -o main

if [ $? -eq 0 ]; then
    echo ---Builded!---
    ./main
    echo -----END------
else
    echo Build Failed!
fi 