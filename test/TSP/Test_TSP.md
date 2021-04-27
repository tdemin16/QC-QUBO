
# First
Test conducted by removing the tabu matrix and mantaining the probability of permutation = 1.

|     x    |     y     |
|----------|-----------|
|7.21362447| 4.72776071|
|7.48093569| 8.15635289|
|9.72339245| 7.86815924|
|2.68741026| 2.43794644|
|6.86692814| 2.86802978|
|8.11962058| 4.72617735|
|9.66730978| 1.6601765 |
|0.48265261| 9.97703727|
|4.85291697| 3.57640542|
|4.64010762| 0.77011666|

## Brute Force solution
Brute Force: [0, 2, 1, 7, 3, 9, 8, 4, 6, 5] 36.2746 \
Calculation time: 0.475904

## D-Wave solution
D-Wave: [2, 1, 4, 3, 9, 0, 6, 5, 7, 8] 49.947659879310976 \
Calculation time: 117.91059374809265

## Hybrid solution
Hybrid: [2, 1, 7, 3, 4, 8, 9, 6, 0, 5] 39.96508710124576 \
Calculation time: 17.37695550918579

## QALS solution
QALS: [9, 4, 5, 2, 6, 3, 1, 0, 8, 7] 53.439\
Calculation time: 38.6439

QALS: [9, 1, 0, 7, 8, 6, 4, 2, 5, 3] 53.619\
Calculation time: 78.8585

QALS: [6, 2, 7, 5, 1, 9, 4, 3, 8, 0] 52.6196\
Calculation time: 88.2038

# Second
Restored tabu search with $\lambda_0 = 3 / 2$ and N = 5.

|     x    |     y     |
|----------|-----------|
|7.65633375| 0.78444412|
|1.74143071| 3.71198309|
|1.49034874| 7.7039062 |
|3.57395371| 5.4498407 |
|2.1242206 | 5.70778437|
|9.25808422| 4.19477601|
|2.99161883| 5.25690292|
|1.34852573| 3.4841625 |
|7.03554598| 3.29164378|
|8.70281554| 9.5195283 |

## Brute Force solution
Brute Force: [0, 5, 9, 2, 4, 7, 1, 6, 3, 8] 30.7255 \
Calculation time: 0.519616

## D-Wave solution
D-Wave: [6, 9, 8, 7, 2, 5, 1, 0, 4, 3] 55.633240896012566 \
Calculation time: 115.58193039894104

## Hybrid solution
Hybrid: [4, 2, 9, 5, 8, 0, 7, 1, 6, 3] 31.256196933482386 \
Calculation time: 16.84376621246338

## QALS solution
QALS: [6, 2, 7, 5, 1, 9, 4, 3, 8, 0] 52.6196 \
Calculation time: 88.2038

QALS: [9, 0, 6, 4, 3, 7, 2, 5, 1, 8] 48.345 \
Calculation time: 57.9059

# Third
Test with n = 5 and same conditions of the second test

|     x    |     y     |
|----------|-----------|
|8.52591888| 3.21378717|
|3.03452934| 5.30414046|
|8.9685915 | 7.96356731|
|4.48992494| 7.32524818|
|4.71858654| 6.57472467|

## Brute Force solution
Brute Force: [0, 1, 4, 3, 2] 18.0643 \
Calculation time: 6.1767e-05

## D-Wave solution
D-Wave: [2, 3, 4, 1, 0] 18.064270687463058 \
Calculation time: 4.9984002113342285

## Hybrid solution
Hybrid: [1, 4, 3, 2, 0] 18.064270688232963 \
Calculation time: 16.947219610214233

## QALS solution
QALS: [4, 0, 1, 2, 3] 22.7656 \
Calculation time: 23.5555


# Fourth
n = 12 which is the maximum dimension for embedding. 
$$
	12^2 = 144 \qquad \text{dimension of QUBO matrix}
$$
$$
	\frac{144^2}{4} = 5184 \qquad \text{number of qubits required for embedding}
$$

|     x    |     y     |
|----------|-----------|
|8.8580261 | 9.19800721|
|3.58473975| 1.65934421|
|4.6050223 | 0.72658109|
|0.01825975| 3.78999345|
|0.16127204| 5.97264417|
|7.03489676| 7.40024074|
|0.86655297| 0.49004868|
|0.35013547| 7.08491828|
|7.02722521| 5.99131315|
|4.32180413| 1.06361234|
|3.81010642| 7.63817388|
|4.45345072| 4.45003767|

## Brute Force solution
Brute Force: [0, 5, 8, 11, 1, 9, 2, 6, 3, 4, 7, 10] 30.5361 \
Calculation time: 62.0376

## D-Wave solution
D-Wave: [0, 8, 2, 7, 11, 9, 5, 6, 10, 4, 3, 1] 68.84281051637674 \
Calculation time: 179.92916774749756

## Hybrid solution
Hybrid: [4, 3, 6, 9, 2, 10, 0, 5, 8, 1, 11, 7] 40.20793498986772
Calculation time: 8.500251531600952

## QALS solution
QALS: [11, 0, 3, 2, 10, 9, 5, 8, 7, 6, 4, 1] 71.5469 \
Calculation time: 59.7533

QALS: [4, 1, 2, 9, 6, 0, 11, 10, 3, 7, 5, 8] 56.0611 \
Calculation time: 115.613


###NUOVI TEST

[[9.94842154 7.98192104]
 [8.60849611 8.00248268]
 [7.7510664  8.92952639]
 [2.48323566 7.71263975]
 [0.79106925 4.7339348 ]
 [5.43068116 8.76787367]
 [2.01200685 7.57374064]
 [8.24493747 4.47630442]
 [6.12997546 7.6483757 ]
 [5.42759893 6.52504111]
 [1.19438987 6.82360516]
 [5.30828759 3.38905631]] 

Hybrid solution
Hybrid:  [ 6  3  8  9  0  1  2  7 11  4 10  5] 35.55181818926448
Calculation time: 14.736042976379395
D-Wave solution
D-Wave: [ 1 10  0  6  9  7 11  4  5  3  8  2] 55.44057099593448
Calculation time: 102.12452530860901
QALS solution
QALS: [7, 11, 6, 0, 8, 4, 10, 5, 9, 3, 1, 2] 50.4034
Calculation time: 114.316
Brute Force solution
Brute Force: [0, 1, 2, 5, 8, 9, 3, 6, 10, 4, 11, 7] 26.2198
Calculation time: 58.8898

#2
Hybrid solution
Hybrid:  [ 6 10  3  5  9  1  8  2  0  7 11  4] 33.35232386675365
Calculation time: 15.023589849472046
D-Wave solution
D-Wave: [ 6  1  5 10 11  7  8  9  4  3  2  0] 52.310869978624545
Calculation time: 159.7624580860138
QALS solution
QALS: [5, 0, 9, 8, 7, 10, 2, 3, 6, 11, 1, 4] 60.3087
Calculation time: 114.775
Brute Force solution
Brute Force: [0, 1, 2, 5, 8, 9, 3, 6, 10, 4, 11, 7] 26.2198
Calculation time: 57.4762

#3
Hybrid solution
Hybrid:  [ 9  8  1  2  0  7 11  5 10  3  6  4] 34.67329723130918
Calculation time: 14.903441905975342
D-Wave solution
D-Wave: [ 4  7  3 11  2  5 10  6  8  9  0  1] 53.38985153698601
Calculation time: 137.01140928268433
QALS solution
QALS: [2, 4, 11, 6, 0, 5, 7, 3, 9, 10, 1, 8] 61.9356
Calculation time: 115.906
Brute Force solution
Brute Force: [0, 1, 2, 5, 8, 9, 3, 6, 10, 4, 11, 7] 26.2198
Calculation time: 55.7427

#4
QALS solution
QALS: [2, 5, 1, 4, 7, 8, 0, 9, 3, 10, 6, 11] 51.1538
Calculation time: 112.663
Brute Force solution
Brute Force: [0, 1, 2, 5, 8, 9, 3, 6, 10, 4, 11, 7] 26.2198
Calculation time: 58.055

#5
QALS solution
QALS: [10, 0, 6, 3, 5, 1, 7, 8, 2, 11, 4, 9] 53.0739
Calculation time: 112.203
Brute Force solution
Brute Force: [0, 1, 2, 5, 8, 9, 3, 6, 10, 4, 11, 7] 26.2198
Calculation time: 56.9255
