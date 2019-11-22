% random ids

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
