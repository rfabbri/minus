% random ids for use with synthdata command

% each line is 0-5116 (number of points in synthdata is 5117)
% each frame is 0-99 (number of frames in synthdata is 100)


% 1) Generate permutation 0-99
% 2) Generate second number avoiding repetition
% 3) Generate third number avoiding repetition

f = [...% OK to have nearby views, since this dataset is completely random
0   36 14
1   33 75
2   77 52
3   55 57
4   92 95
5   30 27
6   73 45
7   17 48
8   45 72
9   31 79
10  28 76
11  90 7
12  68 38
13  52 80
14  76 12
15  21 74
16  53 59
17  46 89
18  43 45
19  27 22
20  65 58
21  82 3
22  13 11
23  19 15
24  23 62
25  49 26
26  50 56
27  66 85
28  63 51
29  22 71
30  75 19
31  0  68
32  34 17
33  15 50
34  96 98
35  70 28
36  37 33       
37  39 42
38  91 70
39  99 40
40  8  35
41  18 13
42  16 87
43  80 69
44  94 32
45  29 60
46  14 77
47  11 61
48  7  2
49  85 67
50  9  24
51  48 21
52  72 99
53  61 92
54  10 47
55  32 53
56  57 34
57  3  18
58  87 6
59  54 31
60  59 40
61  40 64
62  42 30
63  4  90
64  5  94
65  95 4
66  98 9
67  97 83
68  35 55
69  26 44
70  69 65
71  56 82
72  1  78
73  64 10
74  75 73
75  6  86
76  12 66
77  47 46
78  86 88
79  80 26
80  83 5
81  51 63
82  81 1
83  62 49
84  25 98
85  44 97
86  93 36
87  78 54
88  84 41
89  24 81
90  20 8
91  41 23
92  67 0
93  58 16
94  88 93
95  89 20
96  2  37
97  60 91
98  71 29
99  38 84
];


% 1) seq 0-99 shuf shuf
% 2) compare:
%
find(f(:,1) == f(:,3))
find(f(:,1) == f(:,2))
find(f(:,2) == f(:,3))

% remove dup lines
% sort -n | uniq

%
% For each of these indices, 
% keep summing 1 until it returns null
%
% Check for repeated triplets

% for i=1:100
%   while f(i,1) == f(i,2)
%     k = mod(i+1,100);
%     while (k ~= i)
%       if f(k,2) ~= f(i,1)
%         swap f(i,2) and f(k,2)
%         aux = f(i,2);
%         f(i,2) = f(k,2);
%         f(k,2) = aux;
%       end
%     k = mod(i+1,100) + 1;
%     end
%   end
% end


p = [
 0  4035    3020  
 50 1985    3670  
100 1135     820 
150  985    3170   
200  2985   4870    
250   485   2920    
300  2235   3970    
350  2085   4220    
400  3135   1170    
450  3885     20    
500   185     70    
550  1635    320    
600  1535   1970    
650  2685   4720    
700  1585   2120    
750  3185    270    
800  2735      4    
900   635   2220    
950  1835   3470    
1000  785   4370     
1050 2785   2720     
1100  685   4320     
1150  935   1620     
1200  335    120     
1250 1335   2420     
1300 3735   1470     
1350 4535   1220     
1400 4135   2127     
1450 4285    470     
1500 3535   4420     
1550 1185    220     
1600 3385   3920     
1650 1035   4770     
1700 2885   1420     
1750 4335   3870     
1800  235   3570     
1850 3485   4670     
1900 1235   1870     
1950 1785   2070     
2000 4435    520     
2050 4585   1120     
2100 2135   4570     
2150 3935   2520     
2200 3685    620     
2250 1935   4920     
2300 1285   3420     
2350 4835   4970     
2400 1085   3220     
2450 1485   3479     
2500 1435   3770     
2550 2185   1370     
2650 2335   1020     
2700 5085    370     
2750 4385   3320     
2800 4985   2570     
2850 4785   4520     
2900 2435   2370     
2950 3585    770     
3000 2835   2020     
3050  435   4020     
3100 3235   5070     
3150 4885   2970     
3200  735   3120     
3250  835   2870     
3300 4787   1720     
3350   35   1770     
3400 3835   4070     
3450 1435   4620     
3500  135   4470     
3550 5035   2770     
3600   85   1070     
3650 3035    570     
3700 2535   1570     
3750 2035   2270     
3800 2385   3520     
3850 1885    670     
3900 2935   3070     
3950 3635   3270     
4000  535   1920     
4050 1385   1320     
4100 3785   1820     
4200  285   3720     
4250 4635   5020     
4300 4485   1270     
4350 1685   1520     
4400 1735   4820     
4450  585    970     
4500   19   3620     
4550  385   1670     
4600 3085    720     
4650 3985   3370     
4700 1938   3820     
4750 738   2320     
4800 2485   2670     
4850 3335    920     
4900 4085    170     
4950 4235    420     
5000 2585   2170     
5050 4685   2820     
5100 2285   4270     
];
% this is then shuffled, then concat with f, and both shuffled, to generate
% 100-configuration-synthdata
pshuf=[
1750 4335   3870     
4550  385   1670     
1000  785   4370     
2350 4835   4970     
2750 4385   3320     
1650 1035   4770     
3800 2385   3520     
1150  935   1620     
4450  585    970     
4300 4485   1270     
1400 4135   2127     
5000 2585   2170     
4850 3335    920     
2050 4585   1120     
5050 4685   2820     
1200  335    120     
650  2685   4720    
1450 4285    470     
2650 2335   1020     
3850 1885    670     
3600   85   1070     
900   635   2220    
1800  235   3570     
550  1635    320    
3950 3635   3270     
4050 1385   1320     
2850 4785   4520     
700  1585   2120    
4650 3985   3370     
1500 3535   4420     
1850 3485   4670     
750  3185    270    
1250 1335   2420     
3750 2035   2270     
4950 4235    420     
2150 3935   2520     
600  1535   1970    
3400 3835   4070     
3250  835   2870     
250   485   2920    
1100  685   4320     
3100 3235   5070     
1350 4535   1220     
4700 1938   3820     
400  3135   1170    
1900 1235   1870     
1050 2785   2720     
4350 1685   1520     
3900 2935   3070     
2900 2435   2370     
4750 738   2320     
4400 1735   4820     
1600 3385   3920     
4100 3785   1820     
500   185     70    
2400 1085   3220     
4250 4635   5020     
3050  435   4020     
950  1835   3470    
350  2085   4220    
4800 2485   2670     
3000 2835   2020     
200  2985   4870    
3550 5035   2770     
100 1135     820 
4900 4085    170     
5100 2285   4270     
3500  135   4470     
2250 1935   4920     
 0  4035    3020  
2500 1435   3770     
3450 1435   4620     
2300 1285   3420     
3300 4787   1720     
4200  285   3720     
3700 2535   1570     
3150 4885   2970     
300  2235   3970    
2200 3685    620     
800  2735      4    
450  3885     20    
1700 2885   1420     
 50 1985    3670  
4500   19   3620     
1550 1185    220     
2800 4985   2570     
2950 3585    770     
4000  535   1920     
2000 4435    520     
3650 3035    570     
150  985    3170   
1950 1785   2070     
2700 5085    370     
2550 2185   1370     
3200  735   3120     
2100 2135   4570     
1300 3735   1470     
3350   35   1770     
4600 3085    720     
2450 1485   3479     
];


mind=25; 

find(abs(p(:,1) - p(:,3)) < mind)
find(abs(p(:,1) - p(:,2)) < mind)
find(abs(p(:,2) - p(:,3)) < mind)
  

% v=[...
% 0   
% 50  
% 100 
% 150 
% 200 
% 250 
% 300 
% 350 
% 400 
% 450 
% 500 
% 550 
% 600 
% 650 
% 700 
% 750 
% 800 
% 900 
% 950 
% 1000
% 1050
% 1100
% 1150
% 1200
% 1250
% 1300
% 1350
% 1400
% 1450
% 1500
% 1550
% 1600
% 1650
% 1700
% 1750
% 1800
% 1850
% 1900
% 1950
% 2000
% 2050
% 2100
% 2150
% 2200
% 2250
% 2300
% 2350
% 2400
% 2450
% 2500
% 2550
% 2650
% 2700
% 2750
% 2800
% 2850
% 2900
% 2950
% 3000
% 3050
% 3100
% 3150
% 3200
% 3250
% 3300
% 3350
% 3400
% 3450
% 3500
% 3550
% 3600
% 3650
% 3700
% 3750
% 3800
% 3850
% 3900
% 3950
% 4000
% 4050
% 4100
% 4200
% 4250
% 4300
% 4350
% 4400
% 4450
% 4500
% 4550
% 4600
% 4650
% 4700
% 4750
% 4800
% 4850
% 4900
% 4950
% 5000
% 5050
% 5100 ];

% v = mod(v+20,5116);

