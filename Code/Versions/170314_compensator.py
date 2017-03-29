# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 13:44:04 2016

@author: Max Karlsson
"""
from pyteomics import mass
from Printicek import daviPrint
import copy
from math import floor
"""
#DEV ONLY:
from mass_calculator import massCalculator
isolationListPath = 'C:/Users/bjorn.forsstrom/Documents/Max/Wellness project/WellnessIsolationList short dummy.csv'
settings = {'Isolation window offset': 1.25, 'Max fragment charge': 2, 'MS2 tolerance': 7.0, 'MS1 tolerance': 10.0, 'Co-isolation': 1, 'Scheduled run': 0, 'Include b-fragments': 1, 'Fragment peptide calling': 0, 'MS2 quantification': 1, 'Reduce noise': 1, 'Heavy labels': 'K:0C2N, R:0C2N', 'Peaks in isotope clusters': 3, 'Isolation window size': 3.0, 'MS1 quantification': 0, 'IAA': 1}
peptideSettings = massCalculator(isolationListPath, settings)
indexedData = {'TPSAAYLWVGTGASEAEK++': {'MS2': {23.150459: {'y14 l0+': 4964.76123046875, 'y9 l0+': 21054.51953125, 'y13 l1+': 2813.41748046875, 'y11 l1+': 7185.19482421875, 'y12 l0+': 13384.2255859375, 'y13 l0+': 5823.703125, 'y10 l1+': 2935.1923828125, 'y11 l0+': 23552.212890625, 'y10 l0+': 13887.4130859375, 'y7 l0+': 5180.93994140625}, 23.420364: {'b4 h0+': 11028.171875, 'b2 h0+': 5132.03173828125, 'b7 h0+': 11067.611328125, 'b5 h0+': 9854.640625, 'b3 h0+': 5154.6748046875, 'b6 h0+': 11239.2158203125}, 22.632284: {'y9 l0+': 84794.0234375, 'y17 l2++': 11460.509765625, 'y8 l0+': 13237.37890625, 'y14 l1+': 10147.6396484375, 'y12 l0+': 45556.046875, 'y11 l1+': 27692.677734375, 'y10 l1+': 7963.7724609375, 'y16 l0+': 7039.541015625, 'y7 l0+': 13446.1142578125, 'y15 l0+': 6271.62744140625, 'y3 l0+': 6333.60498046875, 'y17 l1++': 26408.1875, 'y17 l0++': 7396.560546875, 'y13 l1+': 12261.728515625, 'y11 l0+': 90954.171875, 'y13 l0+': 42155.796875, 'y10 l0+': 51414.83203125, 'y9 l1+': 11411.6845703125, 'y12 l1+': 18093.998046875, 'y17 l1+': 10729.3046875}, 23.889042: {'b4 h0+': 2865.9638671875, 'b7 h0+': 7352.310546875, 'b5 h0+': 8892.93359375, 'b6 h0+': 5066.56689453125, 'b2 h0+': 6127.92138671875}, 22.915088: {'b4 h0+': 12055.998046875, 'b8 h0+': 8467.876953125, 'b2 h0+': 11455.8271484375, 'b7 h0+': 21357.015625, 'b9 h0+': 12887.126953125, 'b5 h0+': 23561.99609375, 'b3 h0+': 6048.7216796875, 'b6 h0+': 15421.2822265625}, 23.284927: {'y14 l0+': 6012.4990234375, 'y9 l0+': 19583.33984375, 'y17 l1++': 31975.744140625, 'y11 l1+': 5056.68994140625, 'y11 l0+': 20868.37890625, 'y12 l0+': 11190.357421875, 'y10 l0+': 11808.947265625, 'y17 l2++': 13588.5869140625}, 23.195355: {'y14 l0+': 6241.0634765625, 'y9 l0+': 15892.2119140625, 'y11 l1+': 10060.306640625, 'y11 l0+': 21003.5390625, 'y12 l0+': 10591.505859375, 'y13 l0+': 10002.1162109375, 'y10 l0+': 15422.37890625}, 22.517634: {'y14 l0+': 49312.8515625, 'y9 l0+': 205956.5625, 'y5 l0+': 14758.0400390625, 'y14 l1+': 28909.16015625, 'y8 l1+': 6242.47021484375, 'y9 l1+': 35926.0703125, 'y2 l0+': 14592.6103515625, 'y13 l1+': 47572.76953125, 'y17 l1++': 13423.859375, 'y11 l1+': 75801.0234375, 'y11 l0+': 259289.984375, 'y3 l0+': 16466.33984375, 'y10 l0+': 149662.28125, 'y15 l1+': 8733.6591796875, 'y12 l1+': 42775.32421875, 'y12 l0+': 144133.953125, 'y6 l0+': 8143.8134765625, 'y8 l0+': 35570.84375, 'y10 l1+': 32995.4921875, 'y16 l0+': 25947.8984375, 'y4 l0+': 7086.61279296875, 'y15 l0+': 22941.091796875, 'y7 l0+': 38482.734375, 'y17 l0++': 7454.48095703125, 'y13 l0+': 119084.1640625, 'y16 l1+': 13026.826171875}, 23.197613: {'b4 h0+': 7967.2509765625, 'b8 h0+': 7785.8525390625, 'b2 h0+': 6172.7578125, 'b7 h0+': 4899.875, 'b9 h0+': 3009.90869140625, 'b5 h0+': 14862.6494140625, 'b3 h0+': 2668.04345703125, 'b6 h0+': 8473.0244140625}, 22.912815: {'y14 l0+': 12442.6083984375, 'y9 l0+': 44816.03125, 'y8 l0+': 5696.72021484375, 'y14 l1+': 6037.0625, 'y7 l0+': 6571.740234375, 'y12 l0+': 28853.26953125, 'y9 l1+': 6295.20654296875, 'y15 l0+': 5574.16552734375, 'y13 l1+': 8991.599609375, 'y17 l0++': 5505.9228515625, 'y11 l1+': 14385.53515625, 'y11 l0+': 46270.0703125, 'y13 l0+': 17086.8125, 'y10 l0+': 25284.064453125, 'y12 l1+': 7227.36572265625}, 23.922808: {'b8 h0+': 3148.598388671875, 'b2 h0+': 3217.740234375, 'b5 h0+': 8936.5732421875, 'b7 h0+': 3118.074951171875, 'b3 h0+': 5894.84375, 'b6 h0+': 8619.8720703125}, 23.472161: {'b4 h0+': 10823.595703125, 'b6 h0+': 8009.04736328125, 'b3 h0+': 6703.85107421875, 'b7 h0+': 7173.52099609375, 'b5 h0+': 11060.873046875}, 22.962327: {'b4 h0+': 11295.0341796875, 'b8 h0+': 9312.7626953125, 'b2 h0+': 7751.54833984375, 'b7 h0+': 11217.5927734375, 'b9 h0+': 13596.5869140625, 'b5 h0+': 16673.611328125, 'b3 h0+': 8247.220703125, 'b6 h0+': 18409.19140625}, 23.062664: {'y14 l0+': 8159.6533203125, 'y9 l0+': 18917.41796875, 'y8 l0+': 5065.6865234375, 'y14 l1+': 5232.32958984375, 'y12 l0+': 20003.484375, 'y10 l1+': 7220.01953125, 'y13 l1+': 3030.3828125, 'y11 l1+': 11355.0673828125, 'y11 l0+': 33512.81640625, 'y13 l0+': 13826.0419921875, 'y10 l0+': 13046.048828125, 'y12 l1+': 4509.11376953125}, 23.23813: {'y14 l0+': 5392.8955078125, 'y9 l0+': 15908.3515625, 'y17 l1++': 5362.3564453125, 'y11 l1+': 6344.01416015625, 'y11 l0+': 23782.505859375, 'y12 l0+': 10735.9794921875, 'y13 l0+': 9472.232421875, 'y10 l0+': 11835.2392578125}, 22.862694: {'y14 l0+': 11541.8388671875, 'y9 l0+': 40875.80859375, 'y8 l0+': 5375.55029296875, 'y12 l0+': 24935.173828125, 'y9 l1+': 10330.8115234375, 'y7 l0+': 6482.78076171875, 'y13 l1+': 9712.0673828125, 'y11 l1+': 17735.158203125, 'y11 l0+': 40300.0625, 'y13 l0+': 20958.33203125, 'y10 l0+': 28688.8671875, 'y12 l1+': 10139.439453125}, 22.574527: {'b5 h1+': 6159.2978515625, 'b7 h1+': 10774.0693359375, 'b6 h1+': 8703.978515625, 'b5 h0+': 80379.65625, 'b9 h0+': 35024.56640625, 'b7 h0+': 65913.53125, 'b4 h0+': 47684.35546875, 'b8 h0+': 35538.33203125, 'b3 h0+': 13431.662109375, 'b2 h0+': 31927.875, 'b9 h1+': 12908.052734375, 'b8 h1+': 6959.7705078125, 'b10 h0+': 4761.060546875, 'b6 h0+': 63492.921875}, 23.016962: {'y14 l0+': 6572.62158203125, 'y9 l0+': 32253.7578125, 'y13 l1+': 5991.3173828125, 'y11 l1+': 12466.83984375, 'y11 l0+': 37585.4921875, 'y12 l0+': 23724.640625, 'y13 l0+': 18536.603515625, 'y9 l1+': 3027.098876953125, 'y10 l0+': 20308.10546875, 'y7 l0+': 6912.35302734375, 'y12 l1+': 6919.56298828125}, 23.990995: {'y9 l0+': 11763.4892578125, 'y13 l1+': 2819.267822265625, 'y11 l0+': 16853.787109375, 'y12 l0+': 7007.22509765625, 'y13 l0+': 9300.4833984375, 'y10 l1+': 2629.817626953125, 'y10 l0+': 8661.7646484375, 'y12 l1+': 2946.5771484375}, 23.959487: {'b8 h0+': 2778.43994140625, 'b6 h0+': 2697.570556640625, 'b9 h0+': 5105.14013671875, 'b5 h0+': 6353.56494140625}, 23.107711: {'b4 h0+': 5596.6015625, 'b2 h0+': 6000.27001953125, 'b7 h0+': 8224.439453125, 'b9 h0+': 6419.3876953125, 'b5 h0+': 13403.9677734375, 'b6 h0+': 10878.9716796875}, 22.864953: {'b4 h0+': 10941.5888671875, 'b8 h0+': 11116.384765625, 'b2 h0+': 9087.6044921875, 'b7 h0+': 19688.0234375, 'b9 h0+': 10193.02734375, 'b5 h0+': 22120.232421875, 'b3 h0+': 4873.2568359375, 'b6 h0+': 14715.326171875}, 23.064924: {'b4 h0+': 10935.0234375, 'b8 h0+': 4039.529052734375, 'b2 h0+': 7988.7666015625, 'b7 h0+': 13121.5361328125, 'b9 h0+': 6453.25830078125, 'b5 h0+': 11083.1962890625, 'b6 h0+': 13037.7841796875}, 23.152718: {'b4 h0+': 6869.9111328125, 'b8 h0+': 7322.46142578125, 'b2 h0+': 6516.04931640625, 'b7 h0+': 8234.517578125, 'b9 h0+': 5954.47802734375, 'b5 h0+': 13373.0390625, 'b6 h0+': 8087.3994140625}, 23.327706: {'y9 l0+': 20284.537109375, 'y17 l1++': 34975.5625, 'y11 l1+': 6395.61572265625, 'y11 l0+': 22530.146484375, 'y12 l0+': 8449.8994140625, 'y13 l0+': 9690.8466796875, 'y10 l0+': 10947.8740234375, 'y17 l2++': 22776.8828125, 'y12 l1+': 7211.3759765625, 'y17 l1+': 8316.8056640625}, 23.019225: {'b4 h0+': 11961.2138671875, 'b8 h0+': 8788.4873046875, 'b2 h0+': 7538.65869140625, 'b7 h0+': 13907.5263671875, 'b9 h0+': 10372.3486328125, 'b5 h0+': 16113.87109375, 'b6 h0+': 14937.333984375}, 23.957226: {'y14 l0+': 2809.75341796875, 'y12 l0+': 9980.4921875, 'y9 l0+': 13567.1630859375, 'y13 l0+': 9119.26171875, 'y11 l0+': 15385.4248046875, 'y10 l0+': 7233.4287109375}, 23.329964: {'b4 h0+': 6373.9130859375, 'b8 h0+': 5655.77490234375, 'b7 h0+': 5740.88330078125, 'b5 h0+': 10147.6884765625, 'b3 h0+': 6264.6962890625, 'b6 h0+': 11031.962890625}, 23.240388: {'b4 h0+': 5602.08984375, 'b7 h0+': 8834.6865234375, 'b5 h0+': 11086.3310546875, 'b6 h0+': 8658.94921875, 'b2 h0+': 5011.15673828125}, 22.960066: {'y14 l0+': 10134.4189453125, 'y9 l0+': 34771.359375, 'y8 l0+': 7413.35791015625, 'y14 l1+': 3108.365478515625, 'y12 l0+': 26661.80078125, 'y9 l1+': 7776.69482421875, 'y16 l0+': 5923.12255859375, 'y7 l0+': 4948.7265625, 'y13 l1+': 9124.892578125, 'y11 l1+': 14329.333984375, 'y11 l0+': 44540.1796875, 'y13 l0+': 22150.560546875, 'y10 l0+': 20406.873046875, 'y12 l1+': 6211.6572265625}, 22.572267: {'y14 l0+': 36180.88671875, 'y9 l0+': 146540.171875, 'y5 l0+': 12797.978515625, 'y14 l1+': 15991.26953125, 'y9 l1+': 26751.947265625, 'y2 l0+': 7371.0712890625, 'y3 l0+': 16941.630859375, 'y13 l1+': 32177.794921875, 'y17 l1++': 13839.943359375, 'y11 l1+': 57276.640625, 'y11 l0+': 177767.28125, 'y4 l0+': 7410.94287109375, 'y10 l0+': 100166.1953125, 'y15 l1+': 4727.4140625, 'y12 l1+': 34872.67578125, 'y12 l0+': 91897.796875, 'y6 l0+': 5441.56201171875, 'y8 l0+': 20967.11328125, 'y10 l1+': 19234.068359375, 'y16 l0+': 18957.765625, 'y15 l0+': 16724.29296875, 'y7 l0+': 22175.62890625, 'y17 l0++': 8972.064453125, 'y13 l0+': 72237.5, 'y16 l1+': 7879.26513671875}, 23.371637: {'y9 l0+': 19239.091796875, 'y17 l2++': 32082.197265625, 'y11 l0+': 24289.78125, 'y12 l0+': 10206.919921875, 'y13 l0+': 10370.1533203125, 'y10 l0+': 14362.076171875, 'y17 l1++': 60461.6171875, 'y17 l1+': 9182.1396484375}, 22.634544: {'b5 h1+': 5663.4833984375, 'b7 h1+': 6386.38427734375, 'b6 h1+': 6486.1396484375, 'b5 h0+': 56455.984375, 'b9 h0+': 25551.291015625, 'b7 h0+': 40900.96875, 'b4 h0+': 28491.94140625, 'b8 h0+': 21744.76171875, 'b3 h0+': 6724.2216796875, 'b2 h0+': 23700.21875, 'b9 h1+': 5946.521484375, 'b8 h1+': 8667.509765625, 'b6 h0+': 39456.0234375}, 23.287187: {'b4 h0+': 7933.375, 'b2 h0+': 6749.25390625, 'b7 h0+': 9459.8505859375, 'b9 h0+': 6715.56982421875, 'b5 h0+': 14028.0078125, 'b6 h0+': 9750.0732421875}, 23.105453: {'y14 l0+': 5184.09912109375, 'y9 l0+': 22581.68359375, 'y8 l0+': 2804.0234375, 'y12 l0+': 13586.5693359375, 'y9 l1+': 4621.416015625, 'y16 l0+': 3031.02978515625, 'y7 l0+': 4784.71240234375, 'y13 l1+': 6763.78125, 'y11 l1+': 9236.150390625, 'y11 l0+': 32088.423828125, 'y13 l0+': 12844.7548828125, 'y10 l0+': 15707.26953125, 'y12 l1+': 6545.6337890625}, 22.519892: {'b5 h1+': 8510.1044921875, 'b7 h1+': 18370.86328125, 'b6 h1+': 10793.5341796875, 'b5 h0+': 118880.828125, 'b9 h0+': 59258.390625, 'b7 h0+': 104845.0625, 'b4 h0+': 55632.046875, 'b8 h0+': 61762.14453125, 'b3 h0+': 22855.349609375, 'y11 l1+': 7322.544921875, 'b2 h0+': 45708.26171875, 'b9 h1+': 15506.720703125, 'b8 h1+': 16203.345703125, 'b10 h0+': 10781.091796875, 'b6 h0+': 103506.1640625}, 23.993255: {'b4 h0+': 4893.1474609375, 'b9 h0+': 4833.73046875, 'b5 h0+': 8232.1875, 'b6 h0+': 8700.1533203125, 'b2 h0+': 5412.828125}}, 'MS1': {}}, 'SSEDPNEDIVER++': {'MS2': {11.379198: {'y9 l0+': 551251.75, 'y9 l0++': 225360.890625, 'y8 l1+': 347870.46875, 'y8 l0++': 1093021.25, 'y9 l1+': 112770.1640625, 'y10 l1+': 48113.69921875, 'y2 l0+': 156787.390625, 'y3 l0+': 244145.1875, 'y7 l0+': 180803.078125, 'y8 l0+': 2140608.5, 'y5 l0+': 254384.046875, 'y4 l0+': 100115.3671875, 'y6 l0+': 146867.296875, 'y10 l0+': 153993.21875, 'y8 l1++': 153571.3125}, 11.117339: {'y9 l0+': 311915.0, 'y5 l0+': 137951.359375, 'y8 l1+': 216298.59375, 'y5 l1+': 19669.083984375, 'y9 l1+': 76896.8125, 'y9 l0++': 114249.2421875, 'y2 l0+': 89602.5703125, 'y3 l1+': 7290.69482421875, 'y11 l0++': 16931.349609375, 'y3 l0+': 132501.796875, 'y6 l0+': 80285.890625, 'y4 l1+': 9733.4521484375, 'y10 l0+': 76138.6875, 'y8 l0+': 1062837.125, 'y8 l0++': 542676.25, 'y10 l1+': 20789.013671875, 'y9 l1++': 22910.54296875, 'y4 l0+': 77265.4140625, 'y10 l0++': 21695.13671875, 'y7 l1+': 16529.787109375, 'y7 l0+': 112596.7109375, 'y6 l1+': 13464.0322265625, 'y8 l1++': 102424.2578125}, 12.055869: {'y5 l2+': 4293.802734375, 'y6 l0+': 16104.90625, 'y6 l1+': 12115.2119140625, 'y6 l2+': 5125.267578125}, 11.211154: {'y10 l0++': 112025.734375, 'y9 l0+': 1526251.5, 'y8 l0+': 5473128.5, 'y8 l1+': 1068542.875, 'y5 l1+': 74474.7734375, 'y8 l0++': 2696013.25, 'y9 l1+': 323523.5, 'y9 l0++': 572964.125, 'y2 l0+': 488009.59375, 'y9 l1++': 127940.8203125, 'y3 l0+': 616830.4375, 'y7 l1+': 107917.4375, 'y7 l0+': 558316.875, 'y5 l0+': 662936.5625, 'y4 l0+': 396152.21875, 'y8 l1++': 560464.375, 'y6 l0+': 367222.25, 'y10 l0+': 350779.125, 'y10 l1+': 109015.84375}, 12.088127: {'y6 l0+': 6869.287109375}, 11.300961: {'b4 h0+': 1830040.25, 'b9 h0+': 167595.078125, 'b7 h0+': 229160.890625, 'b3 h0+': 610619.75, 'b2 h0+': 1859054.5}, 11.119583: {'b4 h0+': 551865.0625, 'b8 h0+': 24329.0390625, 'b4 h1+': 44379.51171875, 'b2 h0+': 541848.8125, 'b5 h0+': 40053.2109375, 'b10 h0+': 22903.119140625, 'b9 h0+': 48810.7109375, 'b7 h0+': 52218.265625, 'b3 h0+': 146122.203125}, 11.412908: {'y9 l0+': 218619.734375, 'y9 l0++': 79874.25, 'y8 l1+': 150948.03125, 'y8 l0++': 433738.15625, 'y9 l1+': 60482.359375, 'y10 l1+': 21362.267578125, 'y2 l0+': 74596.8203125, 'y9 l1++': 16253.6396484375, 'y3 l0+': 102713.6484375, 'y7 l0+': 85010.3671875, 'y8 l0+': 859961.6875, 'y5 l0+': 130139.8515625, 'y4 l0+': 40676.578125, 'y6 l0+': 61348.11328125, 'y10 l0+': 52391.05859375, 'y8 l1++': 71075.125}, 11.449535: {'y9 l0+': 73113.21875, 'y9 l0++': 24009.59375, 'y8 l1+': 41088.0078125, 'y8 l0++': 123902.71875, 'y9 l1++': 5729.8525390625, 'y9 l1+': 17870.3984375, 'y10 l1+': 5078.1298828125, 'y2 l0+': 14070.1904296875, 'y10 l0++': 4042.301025390625, 'y3 l0+': 29908.85546875, 'y7 l0+': 25400.953125, 'y8 l0+': 260240.515625, 'y5 l0+': 31988.31640625, 'y4 l0+': 16443.572265625, 'y6 l0+': 15670.3837890625, 'y10 l0+': 17220.779296875, 'y8 l1++': 17339.490234375}, 11.074094: {'b4 h0+': 166195.28125, 'b8 h0+': 6907.57421875, 'b4 h1+': 11513.2587890625, 'b2 h0+': 144453.03125, 'b7 h0+': 14693.427734375, 'b10 h0+': 5483.384765625, 'b9 h0+': 12589.818359375, 'b5 h0+': 6380.4892578125, 'b3 h0+': 38007.42578125, 'b6 h0+': 5866.19189453125}, 11.253128: {'y10 l0++': 139242.75, 'y9 l0+': 1804861.25, 'y8 l0+': 5923565.5, 'y9 l1+': 417348.0, 'y8 l1+': 1119235.125, 'y5 l1+': 62852.2265625, 'y8 l0++': 2900906.25, 'y10 l1+': 84322.15625, 'y9 l0++': 639199.0, 'y2 l0+': 452853.25, 'y9 l1++': 126859.9375, 'y3 l0+': 658765.5625, 'y7 l1+': 105424.5234375, 'y7 l0+': 655899.25, 'y11 l0++': 80611.6953125, 'y5 l0+': 686430.125, 'y4 l0+': 360982.5625, 'y6 l0+': 493947.03125, 'y10 l0+': 423678.6875, 'y8 l1++': 440412.4375}, 11.16295: {'y9 l0+': 874120.3125, 'y5 l0+': 367315.625, 'y8 l1+': 524231.75, 'y5 l1+': 57724.15625, 'y9 l1+': 207176.3125, 'y9 l0++': 285075.40625, 'y2 l0+': 227409.015625, 'y3 l0+': 345788.75, 'y11 l0++': 32679.951171875, 'y10 l1++': 29967.9140625, 'y4 l0+': 233300.3125, 'y6 l0+': 259057.078125, 'y10 l0+': 225258.96875, 'y8 l0+': 3171879.75, 'y8 l0++': 1556801.875, 'y8 l1++': 331553.9375, 'y9 l1++': 49397.59765625, 'y7 l0+': 384171.65625, 'y7 l1+': 60178.4765625, 'y10 l1+': 62443.5859375, 'y6 l1+': 39386.203125, 'y10 l0++': 40894.65625}, 11.341432: {'y10 l0++': 89828.375, 'y9 l0+': 1286164.125, 'y8 l0+': 4407595.5, 'y8 l1+': 852765.125, 'y5 l1+': 72621.296875, 'y8 l0++': 1889971.75, 'y9 l1+': 314136.71875, 'y9 l0++': 388042.75, 'y2 l0+': 333922.0, 'y9 l1++': 119503.28125, 'y3 l0+': 557429.25, 'y7 l1+': 70156.4921875, 'y7 l0+': 465385.90625, 'y5 l0+': 535408.5625, 'y4 l0+': 315410.0625, 'y8 l1++': 444176.0625, 'y6 l0+': 312849.6875, 'y10 l0+': 334220.84375, 'y10 l1+': 80225.2578125}, 11.485491: {'b4 h0+': 17331.3984375}, 12.125957: {'y6 l0+': 2647.498779296875}, 11.255374: {'b4 h0+': 1772873.25, 'b7 h0+': 188731.28125, 'b4 h1+': 142172.375, 'b9 h0+': 130897.7578125, 'b2 h0+': 1844218.375, 'b3 h0+': 454737.28125}, 11.451781: {'b4 h0+': 65468.67578125, 'b9 h0+': 6057.63916015625, 'b7 h0+': 6582.69482421875, 'b3 h0+': 14493.212890625, 'b2 h0+': 59696.7578125}, 11.298715: {'y9 l0+': 1495208.25, 'y8 l0+': 4956155.0, 'y8 l1+': 930318.375, 'y8 l0++': 2376688.0, 'y8 l1++': 435021.21875, 'y9 l0++': 578192.8125, 'y4 l0+': 313234.90625, 'y2 l0+': 406848.25, 'y9 l1++': 106463.4375, 'y7 l0+': 537282.875, 'y5 l0+': 701582.625, 'y6 l0+': 384725.53125, 'y3 l0+': 600752.5, 'y10 l0+': 369086.0625, 'y9 l1+': 241981.8125}, 12.238222: {'y5 l2+': 2461.18798828125}, 11.415154: {'b4 h0+': 221661.078125, 'b4 h1+': 15395.2685546875, 'b7 h0+': 17125.544921875, 'b10 h0+': 11580.140625, 'b9 h0+': 21418.255859375, 'b2 h0+': 208930.84375, 'b3 h0+': 59172.32421875}, 11.343677: {'b4 h0+': 1265458.0, 'b9 h0+': 103738.6328125, 'b7 h0+': 76802.3125, 'b3 h0+': 256043.796875, 'b2 h0+': 1140142.5}, 11.07185: {'y10 l0++': 7779.19189453125, 'y9 l0+': 114279.8671875, 'y8 l0+': 384373.28125, 'y8 l1+': 76235.625, 'y5 l1+': 9171.142578125, 'y8 l0++': 173647.6875, 'y9 l1+': 25120.939453125, 'y9 l0++': 36532.265625, 'y2 l0+': 35457.375, 'y9 l1++': 8283.068359375, 'y3 l0+': 40706.44140625, 'y7 l1+': 6975.08056640625, 'y7 l0+': 41705.859375, 'y5 l0+': 44677.16015625, 'y4 l0+': 27941.193359375, 'y8 l1++': 32993.16015625, 'y6 l0+': 30231.5390625, 'y10 l0+': 28970.884765625, 'y10 l1+': 7597.90673828125}, 11.483245: {'y9 l0+': 14483.6337890625, 'y7 l0+': 6021.5869140625, 'y5 l0+': 8805.697265625, 'y8 l0+': 61681.25, 'y8 l1+': 8524.6650390625, 'y6 l0+': 6499.59765625, 'y8 l0++': 27551.541015625, 'y3 l0+': 8133.1064453125, 'y9 l0++': 6951.86572265625, 'y2 l0+': 6699.99267578125}, 11.381444: {'b4 h0+': 615489.75, 'b7 h0+': 57737.52734375, 'b3 h0+': 185798.90625, 'b2 h0+': 680677.25}, 11.165195: {'b4 h0+': 1197730.125, 'b8 h0+': 61509.46875, 'b2 h0+': 1140154.0, 'b7 h0+': 103099.953125, 'b4 h1+': 95105.2734375, 'b9 h0+': 105129.515625, 'b5 h0+': 66890.4140625, 'b3 h0+': 287312.34375, 'b8 h1++': 49301.5, 'b6 h0+': 51303.02734375}, 11.213399: {'b4 h0+': 2016403.25, 'b8 h0+': 105794.109375, 'b7 h0+': 161939.375, 'b4 h1+': 170504.40625, 'b9 h0+': 167960.234375, 'b2 h0+': 2015227.25, 'b3 h0+': 555660.9375}}, 'MS1': {}}, 'SWPAVGNCSSALR++': {'MS2': {}, 'MS1': {}}, 'TVAACNLPIVR++': {'MS2': {11.268874: {'y9 l1++': 5103.326171875}}, 'MS1': {}}, 'FNAVLTNPQGDYDTSTGK++': {'MS2': {18.791707: {'y12 l2+': 8147.03173828125}, 18.318366: {'b4 h0+': 48093.84375, 'b2 h0+': 26928.529296875, 'b7 h0+': 19682.630859375, 'b4 h1+': 6417.17529296875, 'b5 h0+': 23917.0625, 'b3 h0+': 74380.234375, 'b6 h0+': 16948.958984375}, 18.210013: {'b5 h1+': 12718.3515625, 'b7 h1+': 27754.482421875, 'b2 h1+': 5959.50048828125, 'b5 h0+': 110587.140625, 'b4 h1+': 20573.3125, 'b9 h0+': 12486.498046875, 'b7 h0+': 110452.8046875, 'b16 h0+': 7423.009765625, 'b4 h0+': 211545.390625, 'b3 h1+': 19864.837890625, 'b6 h1+': 10887.5439453125, 'b11 h0+': 9223.1767578125, 'b2 h0+': 127055.4921875, 'b10 h0+': 6512.29638671875, 'b3 h0+': 344427.03125, 'b6 h0+': 73014.140625}, 18.147784: {'b5 h1+': 13464.572265625, 'b13 h0+': 9124.4658203125, 'b12 h0+': 5786.0166015625, 'b7 h1+': 16535.955078125, 'b6 h1+': 7304.94140625, 'b8 h0+': 9556.7587890625, 'b5 h0+': 80863.2734375, 'b4 h1+': 13812.90625, 'b9 h0+': 12665.7578125, 'b7 h0+': 84660.265625, 'b16 h0+': 7285.57568359375, 'b4 h0+': 146058.46875, 'b3 h1+': 19924.55859375, 'b11 h0+': 6520.26904296875, 'b2 h0+': 102005.3125, 'b3 h0+': 263760.09375, 'b6 h0+': 57563.52734375}, 18.372982: {'b3 h0+': 8359.158203125}, 18.757932: {'y12 l1++': 9515.8173828125}, 18.261473: {'b5 h1+': 10177.220703125, 'b7 h1+': 11958.177734375, 'b6 h1+': 5255.111328125, 'b8 h0+': 8778.5751953125, 'b5 h0+': 74141.3125, 'b4 h1+': 11947.8125, 'b9 h0+': 10170.708984375, 'b7 h0+': 72125.2421875, 'b14 h0+': 5203.4638671875, 'b4 h0+': 132758.421875, 'b3 h1+': 13435.275390625, 'b11 h0+': 6786.361328125, 'b2 h0+': 83702.015625, 'b10 h0+': 9650.89453125, 'b3 h0+': 229904.859375, 'b6 h0+': 46217.296875}, 18.082757: {'b4 h0+': 60102.4453125, 'b3 h1+': 8483.3251953125, 'b7 h1+': 5952.25537109375, 'b2 h0+': 41776.90234375, 'b7 h0+': 31046.83984375, 'b4 h1+': 6424.603515625, 'b5 h0+': 29249.93359375, 'b3 h0+': 106744.28125, 'b6 h0+': 26541.482421875}}, 'MS1': {}}}
"""
class compensationObject:
    """ An object that contains compensation constants for each fragment.
        Peptide specific.
        
    
    """
    def __init__(self, sequence, charge, includedPeaks, fragments, labels, IAA, peaksInCluster):
        self.sequence = sequence 
        self.charge = charge
        self.includedPeaks = includedPeaks
        #Calculates the number of extra neutrons in label:
        self.__calculate_labelNeutrons(labels)
        #Calculates which isoforms are present in the isolation window:
        self.__calculate_isolatedIsoforms()
        #Calculate MS1 distrbution:
        flatten = lambda l: [item for sublist in l for item in sublist]
        self.MS1clusterPeaks = flatten([('l'+str(i),'h'+str(i)) for i in range(peaksInCluster)])
        self.MS1distribution = self.__calculate_distribution(sequence, self.charge, IAA, 
                                                             includedPeaks=self.MS1clusterPeaks)
        #Calculates which peaks overlap in MS2 fragments:
        self.__calculate_allOverlappingPeaks(fragments, peaksInCluster)
        #Calculates compensation constants:
        self.__calculate_compensationConstants(fragments, IAA)
        
        
        
    def __calculate_composition(self, sequence, charge, IAA):
        """ Returns a dictionary of the atomic composition of a given sequence.
        
        MOVE THIS TO COMPENSATION TOOLS

        """
        #Get composition:
        comp = mass.mass.Composition(sequence=sequence)
        #Add H for each charge:
        comp['H'] += charge
        if IAA:
            cys = sequence.count('C')
            comp['C'] += cys*2
            comp['N'] += cys
            comp['O'] += cys
            comp['H'] += cys*3
        #Make monoisotopic composition:
        monoisotopicComp = mass.mass.Composition({'C[12]': comp['C'],  'C[13]': 0,  'H[1]': comp['H'], 'H[2]': 0, 'N[14]': comp['N'],  'N[15]': 0,  'O[16]':  comp['O'], 'O[17]': 0, 'O[18]': 0, 'S[32]': comp['S'], 'S[33]': 0, 'S[34]': 0})
        return monoisotopicComp
        
    def __calculate_labelNeutrons(self, labels):
        """ Calculates the number of extra neutrons in the K vs the R label.
        
        MOVE THIS TO COMPENSATION TOOLS AND CHANGE NAME TO MORE INTUITIVE
        
        """
        self.KC13label = int(labels.split(',')[0].split(':')[1][0])
        self.RC13label = int(labels.split(' ')[1].split(':')[1][0])
        self.KN15label = int(labels.split(',')[0].split(':')[1][2])
        self.RN15label = int(labels.split(' ')[1].split(':')[1][2])
        self.KlabelNeutrons = self.KC13label + self.KN15label
        self.RlabelNeutrons = self.RC13label + self.RN15label
        self.labelNeutrons =  self.KlabelNeutrons*self.sequence.count('K') + self.RlabelNeutrons*self.sequence.count('R')
        
    def __calculate_isolatedIsoforms(self):
        """ Calculates isoforms that are isolated for both heavy and light
            peaks.
        """
        def compositionFilter(isoforms, span, heavy=False):
            """ Filters compositions to only include them with number of extra
                neutrons in a stated span.
            """
            filteredIsoforms = []
            if heavy:
                minC13 = self.KC13label*self.sequence.count('K')
                minC13 +=self.RC13label*self.sequence.count('R')
                minN15 = self.KN15label*self.sequence.count('K')
                minN15 +=self.RN15label*self.sequence.count('R')
            for composition in isoforms:
                C13 = composition['C[13]']
                H2 = composition['H[2]']
                N15 = composition['N[15]']
                O17 = composition['O[17]']
                O18 = composition['O[18]']
                S33 = composition['S[33]']
                S34 = composition['S[34]']
                extraNeutrons = (C13+H2+N15+O17+O18*2+S33+S34*2)
                if heavy and not (composition['C[13]']>=minC13 and composition['N[15]']>=minN15):
                    continue                    
                if span[0]<=extraNeutrons<=span[1]:
                    filteredIsoforms.append(composition)
            return filteredIsoforms
            
        #Find max extra neutrons:
        try:
            maxNeutrons = max([int(peak[-1]) for peak in self.includedPeaks if peak[0] == 'h'])
        except:
            maxNeutrons = 0
        #Add number of neutrons in the label:
        maxNeutrons += self.labelNeutrons
        
        isoforms = []
        for C13 in range(maxNeutrons+1):
            for H2 in range(maxNeutrons+1):
                for N15 in range(maxNeutrons+1):
                    for O17 in range(maxNeutrons+1):
                        for O18 in range(floor(maxNeutrons/2)+1):
                            for S33 in range(maxNeutrons+1):
                                for S34 in range(floor(maxNeutrons/2)+1):
                                    #Peak number:
                                    extraNeutrons = (C13+H2+N15+O17+O18*2+S33+S34*2)
                                    #If the total number of extra neutrons 
                                    #exceed the maximum number of extra 
                                    #neutrons (peak numbers) that are included
                                    #in the isolation window, ignore:
                                    if maxNeutrons<extraNeutrons or extraNeutrons<0: 
                                        continue
                                    #Make relative composition for this 
                                    #specific composition of atoms:
                                    relComp = {}
                                    relComp['C[12]'] = -C13
                                    relComp['C[13]'] = C13
                                    relComp['H[1]'] = -H2
                                    relComp['H[2]'] = H2
                                    relComp['N[14]'] = -N15
                                    relComp['N[15]'] = N15
                                    relComp['O[16]'] = -(O17+O18)
                                    relComp['O[17]'] = O17
                                    relComp['O[18]'] = O18
                                    relComp['S[32]'] = -(S33+S34)
                                    relComp['S[33]'] = S33
                                    relComp['S[34]'] = S34
                                    
                                    isoforms.append(relComp)
        #Get spans of allowed number of extra neutrons for heavy and light sets
        lightNeutronSpan = (0, max([int(peak[-1]) for peak in self.includedPeaks if peak[0] == 'l']))
        heavyNeutronSpan = (self.labelNeutrons, maxNeutrons)
        #Filter isoforms to light and heavy:
        self.lightIsoforms = compositionFilter(isoforms, lightNeutronSpan)
        self.heavyIsoforms = compositionFilter(isoforms, heavyNeutronSpan, heavy=True)
        
    def __calculate_overlappingPeak(self, fragment):
        """ Calculates which peak the inputted peak overlaps with. Makes the
            assumption that peaks with the same number of extra neutrons are
            not resolved. It also assumes that the light monoisotopic peak is 
            included in the isolation window.
            Overlapping peak compensation is only done for peaks that are 
            included in the includedPeaks variable, that is peaks that are 
            within the isolation window.
        """
        fragmentNumber, sequence, charge, heavyLightIdentifyer = self.__getFragmentInfo(fragment)
        peakNumber = int(fragment.strip('+')[-1])
        #Count the number of K and R in sequence:
        numberK = sequence.count('K')
        numberR = sequence.count('R')
        neutronsInLabel = numberK*self.KlabelNeutrons+numberR*self.RlabelNeutrons
        #Pair with potential overlapping peak:
        if heavyLightIdentifyer == 'l':
            overlappingPeak = 'h' + str(peakNumber - neutronsInLabel)
        elif heavyLightIdentifyer == 'h':
            overlappingPeak = 'l' + str(peakNumber + neutronsInLabel)
        if overlappingPeak in self.includedPeaks:
            overlappingPeak = 'y'+str(fragmentNumber)+' '+overlappingPeak+'+'*charge
            return overlappingPeak

        return None
    
    def __calculate_overlappingPrecursorPeak(self, fragment):
        """ Calculates which peak the inputted peak overlaps with. Makes the
            assumption that peaks with the same number of extra neutrons are
            not resolved. It also assumes that the light monoisotopic peak is 
            included in the isolation window.
            Overlapping peak compensation is only done for peaks that are 
            included in the includedPeaks variable, that is peaks that are 
            within the isolation window.
        """
        peakNumber = int(fragment.strip('+')[-1])
        heavyLightIdentifyer = fragment[0]
        #Count the number of K and R in sequence:
        numberK = self.sequence.count('K')
        numberR = self.sequence.count('R')
        neutronsInLabel = numberK*self.KlabelNeutrons+numberR*self.RlabelNeutrons
        #Pair with potential overlapping peak:
        if heavyLightIdentifyer == 'l':
            overlappingPeak = 'h' + str(peakNumber - neutronsInLabel)
        elif heavyLightIdentifyer == 'h':
            overlappingPeak = 'l' + str(peakNumber + neutronsInLabel)
        
        if overlappingPeak in self.MS1clusterPeaks and not int(overlappingPeak[1:])<0 :
            return overlappingPeak
            
        return None
        
        
    def __calculate_allOverlappingPeaks(self, fragments, peaksInCluster):
        """ Calculates peak overlaps for each and every fragment. b-fragments 
            are ignored as they always overlap and should be ignored in 
            co-isolated runs.
        """
        
        self.overlappingPeaks = {}
        #MS2 peaks:
        for fragment in fragments:
            fragmentType = fragment[0]
            #Skip b-fragments:
            if fragmentType == 'b':
                continue
            elif fragmentType == 'y':
                self.overlappingPeaks[fragment] = self.__calculate_overlappingPeak(fragment)
        
        #MS1 peaks:
        #Function for flattening two-dimensional list:
        flatten = lambda l: [item for sublist in l for item in sublist]
        for fragment in flatten([('l'+str(i),'h'+str(i)) for i in range(peaksInCluster)]):
            self.overlappingPeaks[fragment] = self.__calculate_overlappingPrecursorPeak(fragment)

    def __getFragmentInfo(self,fragment):
            """ Gets information from a fragment identifyer string.
            """
            fragmentNumber = int(fragment.split(' ')[0][1:])
            sequence = self.sequence.strip('+')[len(self.sequence.strip('+'))-fragmentNumber:]
            #Count charge
            charge = fragment.count('+')
            #Find identifyer 'h' or 'l':
            heavyLightIdentifyer = fragment.split(' ')[1][0]
            return fragmentNumber, sequence, charge, heavyLightIdentifyer
            
    def __calculate_compensationConstants(self, fragments, IAA):
        """ Calculates compensation constants that can be used to multiply with
            a called MS2 peak to get its MS1 distribution pattern.
        """
        
            
        #Take away '+' if present in sequence:
        precursorSequence = self.sequence.strip('+')
        #Calculate precursor distribution:
        self.precursorDistribution = self.__calculate_distribution(precursorSequence, self.charge, IAA)
        self.heavyPrecursorDistribution = self.__calculate_distribution(precursorSequence, self.charge, IAA, heavy=True)
        
        self.compensationConstants = {}
        
        fragments = [fragment.split(' ')[0]+' '+fragment.split(' ')[1][0]+'+'*fragment.count('+') for fragment in fragments]
        for fragment in sorted(fragments):
            if fragment in self.compensationConstants.keys():
                continue
            fragmentType = fragment[0]
            #Skip b-fragments:
            if fragmentType == 'b':
                continue
            elif fragmentType == 'y':
                fragmentNumber, fragmentSequence, charge, heavyLightIdentifyer = self.__getFragmentInfo(fragment)
                #Get the proper reference distribution set:
                if heavyLightIdentifyer == 'h':
                    heavy=True
                    refPrecursorDistribution = self.heavyPrecursorDistribution
                elif heavyLightIdentifyer == 'l':
                    heavy=False
                    refPrecursorDistribution = self.precursorDistribution
                
                #Calculate distribution:
                distribution = self.__calculate_distribution(fragmentSequence, charge, IAA, heavy=heavy)
                #Calculate compensation constant:
                for peakNumber, abundance in enumerate(distribution):
                    compensationConstant = refPrecursorDistribution[peakNumber]/abundance
                    identifyer = fragment.split(' ')[0] + ' ' + heavyLightIdentifyer+str(peakNumber)+'+'*fragment.count('+')
                    self.compensationConstants[identifyer] = compensationConstant
        
    def __calculate_distribution(self, sequence, charge, IAA, heavy=False, includedPeaks=[]):
        """ Calculates the peak abundance distribution for a specified sequence
            and charge. Can currently only be used for heavy labeled sequences
            where C and N are isotope labeled.
        """
        #If no custom set of included peaks are set, use isolated peaks:
        if len(includedPeaks)==0:
            includedPeaks = self.includedPeaks
        #Distribution among the peaks included in the isolation window. +1 due
        #to that the monoisotopic peak has no extra neutrons:
        if heavy:
            identifyer = 'h'
            isoforms = self.heavyIsoforms
            minC13 = self.KC13label*self.sequence.count('K')
            minC13 +=self.RC13label*self.sequence.count('R')
            minN15 = self.KN15label*self.sequence.count('K')
            minN15 +=self.RN15label*self.sequence.count('R')
        else:
            identifyer = 'l'
            isoforms = self.lightIsoforms
            minC13 = 0
            minN15 = 0
        numberIsolatedPeaks = sum(1 for peak in includedPeaks if peak[0] == identifyer)
        distribution = [0]*numberIsolatedPeaks
        #Calculate composition of ion:
        composition = self.__calculate_composition(sequence, charge, IAA)
        composition['H[1]']+=charge
        #Parse through every possible isoform that contributes to the selected
        #isotope cluster:
        for relativeIsoform in isoforms:
            tempComp = copy.deepcopy(composition)
            extraNeutrons = 0
            #Protect if miscleavage. If miscleaved K or R is removed in 
            #fragmentation, minimum N and C is changed:
            if not sequence.count('K') == self.sequence.count('K') or not sequence.count('R') == self.sequence.count('R'):
                lostK = self.sequence.count('K')-sequence.count('K')
                lostR = self.sequence.count('R')-sequence.count('R')
                tempComp['C[12]'] += (lostK*self.KC13label+lostR*self.RC13label)
                tempComp['C[13]'] -= (lostK*self.KC13label+lostR*self.RC13label)
                tempComp['N[14]'] += (lostK*self.KN15label+lostR*self.RN15label)
                tempComp['N[15]'] -= (lostK*self.KN15label+lostR*self.RN15label)
                #Compensate peak index for lost labels:
                extraNeutrons -= (lostK*self.KC13label+lostR*self.RC13label+lostK*self.KN15label+lostR*self.RN15label)
            tempComp['C[12]'] -= relativeIsoform['C[13]']
            tempComp['C[13]'] += relativeIsoform['C[13]']
            tempComp['H[1]'] -= relativeIsoform['H[2]']
            tempComp['H[2]'] += relativeIsoform['H[2]']
            tempComp['N[14]'] -= relativeIsoform['N[15]']
            tempComp['N[15]'] += relativeIsoform['N[15]']
            tempComp['O[16]'] -= (relativeIsoform['O[17]']+relativeIsoform['O[18]'])
            tempComp['O[17]'] += relativeIsoform['O[17]']
            tempComp['O[18]'] += relativeIsoform['O[18]']
            tempComp['S[32]'] -= (relativeIsoform['S[33]']+relativeIsoform['S[34]'])
            tempComp['S[33]'] += relativeIsoform['S[33]']
            tempComp['S[34]'] += relativeIsoform['S[34]']
            extraNeutrons += relativeIsoform['C[13]'] + relativeIsoform['H[2]'] + relativeIsoform['N[15]'] + relativeIsoform['O[17]'] + relativeIsoform['O[18]']*2 + relativeIsoform['S[33]'] + relativeIsoform['S[34]']*2
            extraNeutrons -= (minC13+minN15)
            
            #Continue if any element is under zero:
            if any([element<0 for element in tempComp.values()]): 
                continue
            #Add abundance of the specific composition:
            distribution[extraNeutrons] += mass.mass.isotopic_composition_abundance(tempComp)
        #Normalisation of distribution to the abundances that were isolated:
        distSum = sum(distribution)
        distribution = [peakAbundance/distSum for peakAbundance in distribution]
        if distribution[0] == 0:
            daviPrint('Abundance zero error.', pause=True)
        distribution = [peak for peak in distribution if peak>0]
        return distribution

def create_compensationLibrary(peptideSettings, settings):
    """ Creates a library of compensation objects that can be used to 
        compensate intensities.
    """
    compensationLibrary = {}
    daviPrint('Calculating compensation constants...', line=True)
    #Parse through peptides:
    for peptide in peptideSettings.values():
        #Get peptide specifics:
        sequence = peptide.sequence
        charge = peptide.charge
        includedPeaks = peptide.includedPeaks
        fragments = peptide.fragments.keys()
        #Print sequence:
        daviPrint(sequence)
        #Create compensation object:
        compObj = compensationObject(sequence, charge, includedPeaks, fragments, settings['Heavy labels'], settings['IAA'], settings['Peaks in isotope clusters'])
        compensationLibrary[sequence+'+'*charge] = compObj
        
    return compensationLibrary

def findMaxMS2Intensity(spectrum, peak):
    """ Finds the peak with the most intensity in the light cluster in the 
        spectrum.
    """
    maxPeak = None
    maxIntensity = 0
    #Get the identifyer of the cluster (e.g. 'y5'):
    clusterID = peak.split(' ')[0]
    charge = peak.count('+')
    for refPeak, intensity in zip(spectrum, spectrum.values()):
        if not peak.startswith(clusterID):
            continue
        heavyLightIdentifyer = refPeak.split(' ')[1][0]
        if heavyLightIdentifyer=='h':
            continue
        if not charge == refPeak.count('+'):
            continue
        if intensity>maxIntensity:
            maxPeak = peak
            maxIntensity = intensity
    return maxPeak, maxIntensity

def calculateTheoreticalMS2Intensity(peak, refPeak, refIntensity, precursorDistribution):
    """ Calculated the theoretical intensity from the light peak with the 
        highest intensity.
    """
    #Get peak numbers to find right abundance:
    peakNumber = int(peak.strip('+')[-1])
    refPeakNumber = int(refPeak.strip('+')[-1])
    peakAbundance = precursorDistribution[peakNumber]
    refPeakAbundance = precursorDistribution[refPeakNumber]
    theoreticalIntensity = peakAbundance*(refIntensity/refPeakAbundance)
    return theoreticalIntensity
    
def compensateMS2data(MS2data, peptide, compObj, compensatedData, compensationConstants):
    for RT in MS2data:
            spectrum = MS2data[RT]
            for peak in spectrum:
                #Ignore if b-fragment:
                if peak[0] == 'b':
                    continue
                overlappingPeak = compObj.overlappingPeaks[peak]
                #Check if there is a overlapping peak:
                if not overlappingPeak==None:
                    heavyLightIdentifyer = peak.split(' ')[1][0]
                    if heavyLightIdentifyer == 'l':
                        #Rename peak to the overlapping peak. Always prioritize
                        #heavy peaks as these should constitute a larger 
                        #intensity:
                        peak.replace('l','h')
                        heavyLightIdentifyer = 'h'
                        
                    if heavyLightIdentifyer == 'h':
                        #Find largest light peak:
                        maxPeak, maxIntensity = findMaxMS2Intensity(spectrum, peak)
                        #If no light peaks were found continue:
                        if maxPeak == None:
                            continue
                        #Calculate theoretical intensity of overlapping peak:
                        theoreticalIntensity = calculateTheoreticalMS2Intensity(
                                overlappingPeak, maxPeak, maxIntensity, 
                                compObj.precursorDistribution)
                        #Subtract theoretical overlapping intensity from peak
                        #intensity:
                        intensity = spectrum[peak]-theoreticalIntensity
                else:
                    intensity = spectrum[peak]
                
                #Multiply peak intensity with compensation constant:
                intensity *= compensationConstants[peak]
                #Check if peak intensity is less than zero. If true, ignore:
                if intensity<=0:
                    continue
                
                #Add compensated peak to new data variable:
                try:
                    compensatedData[peptide]['MS2'][RT][peak] = intensity
                except:
                    compensatedData[peptide]['MS2'][RT] = {}
                    compensatedData[peptide]['MS2'][RT][peak] = intensity
    
    return compensatedData

def findMaxMS1Intensity(spectrum, peak):
    """ Finds the peak with the most intensity in the light cluster in the 
        spectrum.
    """
    maxPeak = None
    maxIntensity = 0
    for refPeak, intensity in zip(spectrum, spectrum.values()):
        heavyLightIdentifyer = refPeak[0]
        if heavyLightIdentifyer=='h':
            continue
        if intensity>maxIntensity:
            maxPeak = peak
            maxIntensity = intensity
    return maxPeak, maxIntensity

def calculateTheoreticalMS1Intensity(peak, refPeak, refIntensity, precursorDistribution):
    """ Calculated the theoretical intensity from the light peak with the 
        highest intensity.
    """
    #Get peak numbers to find right abundance:
    peakNumber = int(peak.strip('+')[-1])
    refPeakNumber = int(refPeak.strip('+')[-1])
    print(peakNumber, precursorDistribution, peak, refPeak)
    peakAbundance = precursorDistribution[peakNumber]
    refPeakAbundance = precursorDistribution[refPeakNumber]
    theoreticalIntensity = peakAbundance*(refIntensity/refPeakAbundance)
    return theoreticalIntensity

def compensateMS1data(MS1data, peptide, compObj, compensatedData):
    for RT in MS1data:
            spectrum = MS1data[RT]
            print(spectrum)
            for peak in spectrum:
                overlappingPeak = compObj.overlappingPeaks[peak]
                #Check if there is a overlapping peak:
                if not overlappingPeak==None:
                    heavyLightIdentifyer = peak[0]
                    if heavyLightIdentifyer == 'l':
                        #Rename peak to the overlapping peak. Always prioritize
                        #heavy peaks as these should constitute a larger 
                        #intensity:
                        peak.replace('l','h')
                        heavyLightIdentifyer = 'h'
                        
                    if heavyLightIdentifyer == 'h':
                        #Find largest light peak:
                        maxPeak, maxIntensity = findMaxMS1Intensity(spectrum, peak)
                        #If no light peaks were found continue:
                        if maxPeak == None:
                            continue
                        #Calculate theoretical intensity of overlapping peak:
                        theoreticalIntensity = calculateTheoreticalMS1Intensity(
                                overlappingPeak, maxPeak, maxIntensity, 
                                compObj.MS1distribution)
                        #Subtract theoretical overlapping intensity from peak
                        #intensity:
                        intensity = spectrum[peak]-theoreticalIntensity
                else:
                    intensity = spectrum[peak]
                
                #Check if peak intensity is less than zero. If true, ignore:
                if intensity<=0:
                    continue
                
                #Add compensated peak to new data variable:
                try:
                    compensatedData[peptide]['MS1'][RT][peak] = intensity
                except:
                    compensatedData[peptide]['MS1'][RT] = {}
                    compensatedData[peptide]['MS1'][RT][peak] = intensity
    
    return compensatedData


def compensateData(compensationLibrary, indexedData, settings):
    compensatedData = {}
    #Create compensation library:    
    for peptide in indexedData:
        if peptide == 'TIC':
            continue
        #Add MS1 data to compensated data:
        MS1data = indexedData[peptide]['MS1']
        MS2data = indexedData[peptide]['MS2']
        
        compensatedData[peptide] = {'MS1':{}, 'MS2':{}}
        #Get compensation object:
        compObj = compensationLibrary[peptide]
        compensationConstants = compObj.compensationConstants
        
        if settings['MS2 compensation']:
            #Compensate MS2data
            compensatedData = compensateMS2data(MS2data, peptide, compObj, compensatedData, compensationConstants)
        
        if settings['MS1 compensation']:
            #Compensate MS1data
            compensatedData = compensateMS1data(MS1data, peptide, compObj, compensatedData)
        
        
        #MAKE MS1 COMPENSATION-------------------------------------------------------------------------------------
                    
    return compensatedData
                    
            
def placeHolder():
    """ 
    >>> from mass_calculator import massCalculator
    >>> isolationListPath = 'C:/Users/max.karlsson/Documents/Max/Wellness project/WellnessIsolationList short dummy.csv'
    >>> settings = {'Isolation window offset': 1.25, 'Max fragment charge': 2, 'MS2 tolerance': 7.0, 'MS1 tolerance': 10.0, 'Co-isolation': 1, 'Scheduled run': 0, 'Include b-fragments': 1, 'Fragment peptide calling': 0, 'MS2 quantification': 1, 'Reduce noise': 1, 'Heavy labels': 'K:0C2N, R:0C2N', 'Peaks in isotope clusters': 3, 'Isolation window size': 3.0, 'MS1 quantification': 0, 'IAA': 1}
    >>> peptideSettings = massCalculator(isolationListPath, settings)
    >>> indexedData = {'TPSAAYLWVGTGASEAEK++': {'MS2': {23.150459: {'y14 l0+': 4964.76123046875, 'y9 l0+': 21054.51953125, 'y13 l1+': 2813.41748046875, 'y11 l1+': 7185.19482421875, 'y12 l0+': 13384.2255859375, 'y13 l0+': 5823.703125, 'y10 l1+': 2935.1923828125, 'y11 l0+': 23552.212890625, 'y10 l0+': 13887.4130859375, 'y7 l0+': 5180.93994140625}, 23.420364: {'b4 h0+': 11028.171875, 'b2 h0+': 5132.03173828125, 'b7 h0+': 11067.611328125, 'b5 h0+': 9854.640625, 'b3 h0+': 5154.6748046875, 'b6 h0+': 11239.2158203125}, 22.632284: {'y9 l0+': 84794.0234375, 'y17 l2++': 11460.509765625, 'y8 l0+': 13237.37890625, 'y14 l1+': 10147.6396484375, 'y12 l0+': 45556.046875, 'y11 l1+': 27692.677734375, 'y10 l1+': 7963.7724609375, 'y16 l0+': 7039.541015625, 'y7 l0+': 13446.1142578125, 'y15 l0+': 6271.62744140625, 'y3 l0+': 6333.60498046875, 'y17 l1++': 26408.1875, 'y17 l0++': 7396.560546875, 'y13 l1+': 12261.728515625, 'y11 l0+': 90954.171875, 'y13 l0+': 42155.796875, 'y10 l0+': 51414.83203125, 'y9 l1+': 11411.6845703125, 'y12 l1+': 18093.998046875, 'y17 l1+': 10729.3046875}, 23.889042: {'b4 h0+': 2865.9638671875, 'b7 h0+': 7352.310546875, 'b5 h0+': 8892.93359375, 'b6 h0+': 5066.56689453125, 'b2 h0+': 6127.92138671875}, 22.915088: {'b4 h0+': 12055.998046875, 'b8 h0+': 8467.876953125, 'b2 h0+': 11455.8271484375, 'b7 h0+': 21357.015625, 'b9 h0+': 12887.126953125, 'b5 h0+': 23561.99609375, 'b3 h0+': 6048.7216796875, 'b6 h0+': 15421.2822265625}, 23.284927: {'y14 l0+': 6012.4990234375, 'y9 l0+': 19583.33984375, 'y17 l1++': 31975.744140625, 'y11 l1+': 5056.68994140625, 'y11 l0+': 20868.37890625, 'y12 l0+': 11190.357421875, 'y10 l0+': 11808.947265625, 'y17 l2++': 13588.5869140625}, 23.195355: {'y14 l0+': 6241.0634765625, 'y9 l0+': 15892.2119140625, 'y11 l1+': 10060.306640625, 'y11 l0+': 21003.5390625, 'y12 l0+': 10591.505859375, 'y13 l0+': 10002.1162109375, 'y10 l0+': 15422.37890625}, 22.517634: {'y14 l0+': 49312.8515625, 'y9 l0+': 205956.5625, 'y5 l0+': 14758.0400390625, 'y14 l1+': 28909.16015625, 'y8 l1+': 6242.47021484375, 'y9 l1+': 35926.0703125, 'y2 l0+': 14592.6103515625, 'y13 l1+': 47572.76953125, 'y17 l1++': 13423.859375, 'y11 l1+': 75801.0234375, 'y11 l0+': 259289.984375, 'y3 l0+': 16466.33984375, 'y10 l0+': 149662.28125, 'y15 l1+': 8733.6591796875, 'y12 l1+': 42775.32421875, 'y12 l0+': 144133.953125, 'y6 l0+': 8143.8134765625, 'y8 l0+': 35570.84375, 'y10 l1+': 32995.4921875, 'y16 l0+': 25947.8984375, 'y4 l0+': 7086.61279296875, 'y15 l0+': 22941.091796875, 'y7 l0+': 38482.734375, 'y17 l0++': 7454.48095703125, 'y13 l0+': 119084.1640625, 'y16 l1+': 13026.826171875}, 23.197613: {'b4 h0+': 7967.2509765625, 'b8 h0+': 7785.8525390625, 'b2 h0+': 6172.7578125, 'b7 h0+': 4899.875, 'b9 h0+': 3009.90869140625, 'b5 h0+': 14862.6494140625, 'b3 h0+': 2668.04345703125, 'b6 h0+': 8473.0244140625}, 22.912815: {'y14 l0+': 12442.6083984375, 'y9 l0+': 44816.03125, 'y8 l0+': 5696.72021484375, 'y14 l1+': 6037.0625, 'y7 l0+': 6571.740234375, 'y12 l0+': 28853.26953125, 'y9 l1+': 6295.20654296875, 'y15 l0+': 5574.16552734375, 'y13 l1+': 8991.599609375, 'y17 l0++': 5505.9228515625, 'y11 l1+': 14385.53515625, 'y11 l0+': 46270.0703125, 'y13 l0+': 17086.8125, 'y10 l0+': 25284.064453125, 'y12 l1+': 7227.36572265625}, 23.922808: {'b8 h0+': 3148.598388671875, 'b2 h0+': 3217.740234375, 'b5 h0+': 8936.5732421875, 'b7 h0+': 3118.074951171875, 'b3 h0+': 5894.84375, 'b6 h0+': 8619.8720703125}, 23.472161: {'b4 h0+': 10823.595703125, 'b6 h0+': 8009.04736328125, 'b3 h0+': 6703.85107421875, 'b7 h0+': 7173.52099609375, 'b5 h0+': 11060.873046875}, 22.962327: {'b4 h0+': 11295.0341796875, 'b8 h0+': 9312.7626953125, 'b2 h0+': 7751.54833984375, 'b7 h0+': 11217.5927734375, 'b9 h0+': 13596.5869140625, 'b5 h0+': 16673.611328125, 'b3 h0+': 8247.220703125, 'b6 h0+': 18409.19140625}, 23.062664: {'y14 l0+': 8159.6533203125, 'y9 l0+': 18917.41796875, 'y8 l0+': 5065.6865234375, 'y14 l1+': 5232.32958984375, 'y12 l0+': 20003.484375, 'y10 l1+': 7220.01953125, 'y13 l1+': 3030.3828125, 'y11 l1+': 11355.0673828125, 'y11 l0+': 33512.81640625, 'y13 l0+': 13826.0419921875, 'y10 l0+': 13046.048828125, 'y12 l1+': 4509.11376953125}, 23.23813: {'y14 l0+': 5392.8955078125, 'y9 l0+': 15908.3515625, 'y17 l1++': 5362.3564453125, 'y11 l1+': 6344.01416015625, 'y11 l0+': 23782.505859375, 'y12 l0+': 10735.9794921875, 'y13 l0+': 9472.232421875, 'y10 l0+': 11835.2392578125}, 22.862694: {'y14 l0+': 11541.8388671875, 'y9 l0+': 40875.80859375, 'y8 l0+': 5375.55029296875, 'y12 l0+': 24935.173828125, 'y9 l1+': 10330.8115234375, 'y7 l0+': 6482.78076171875, 'y13 l1+': 9712.0673828125, 'y11 l1+': 17735.158203125, 'y11 l0+': 40300.0625, 'y13 l0+': 20958.33203125, 'y10 l0+': 28688.8671875, 'y12 l1+': 10139.439453125}, 22.574527: {'b5 h1+': 6159.2978515625, 'b7 h1+': 10774.0693359375, 'b6 h1+': 8703.978515625, 'b5 h0+': 80379.65625, 'b9 h0+': 35024.56640625, 'b7 h0+': 65913.53125, 'b4 h0+': 47684.35546875, 'b8 h0+': 35538.33203125, 'b3 h0+': 13431.662109375, 'b2 h0+': 31927.875, 'b9 h1+': 12908.052734375, 'b8 h1+': 6959.7705078125, 'b10 h0+': 4761.060546875, 'b6 h0+': 63492.921875}, 23.016962: {'y14 l0+': 6572.62158203125, 'y9 l0+': 32253.7578125, 'y13 l1+': 5991.3173828125, 'y11 l1+': 12466.83984375, 'y11 l0+': 37585.4921875, 'y12 l0+': 23724.640625, 'y13 l0+': 18536.603515625, 'y9 l1+': 3027.098876953125, 'y10 l0+': 20308.10546875, 'y7 l0+': 6912.35302734375, 'y12 l1+': 6919.56298828125}, 23.990995: {'y9 l0+': 11763.4892578125, 'y13 l1+': 2819.267822265625, 'y11 l0+': 16853.787109375, 'y12 l0+': 7007.22509765625, 'y13 l0+': 9300.4833984375, 'y10 l1+': 2629.817626953125, 'y10 l0+': 8661.7646484375, 'y12 l1+': 2946.5771484375}, 23.959487: {'b8 h0+': 2778.43994140625, 'b6 h0+': 2697.570556640625, 'b9 h0+': 5105.14013671875, 'b5 h0+': 6353.56494140625}, 23.107711: {'b4 h0+': 5596.6015625, 'b2 h0+': 6000.27001953125, 'b7 h0+': 8224.439453125, 'b9 h0+': 6419.3876953125, 'b5 h0+': 13403.9677734375, 'b6 h0+': 10878.9716796875}, 22.864953: {'b4 h0+': 10941.5888671875, 'b8 h0+': 11116.384765625, 'b2 h0+': 9087.6044921875, 'b7 h0+': 19688.0234375, 'b9 h0+': 10193.02734375, 'b5 h0+': 22120.232421875, 'b3 h0+': 4873.2568359375, 'b6 h0+': 14715.326171875}, 23.064924: {'b4 h0+': 10935.0234375, 'b8 h0+': 4039.529052734375, 'b2 h0+': 7988.7666015625, 'b7 h0+': 13121.5361328125, 'b9 h0+': 6453.25830078125, 'b5 h0+': 11083.1962890625, 'b6 h0+': 13037.7841796875}, 23.152718: {'b4 h0+': 6869.9111328125, 'b8 h0+': 7322.46142578125, 'b2 h0+': 6516.04931640625, 'b7 h0+': 8234.517578125, 'b9 h0+': 5954.47802734375, 'b5 h0+': 13373.0390625, 'b6 h0+': 8087.3994140625}, 23.327706: {'y9 l0+': 20284.537109375, 'y17 l1++': 34975.5625, 'y11 l1+': 6395.61572265625, 'y11 l0+': 22530.146484375, 'y12 l0+': 8449.8994140625, 'y13 l0+': 9690.8466796875, 'y10 l0+': 10947.8740234375, 'y17 l2++': 22776.8828125, 'y12 l1+': 7211.3759765625, 'y17 l1+': 8316.8056640625}, 23.019225: {'b4 h0+': 11961.2138671875, 'b8 h0+': 8788.4873046875, 'b2 h0+': 7538.65869140625, 'b7 h0+': 13907.5263671875, 'b9 h0+': 10372.3486328125, 'b5 h0+': 16113.87109375, 'b6 h0+': 14937.333984375}, 23.957226: {'y14 l0+': 2809.75341796875, 'y12 l0+': 9980.4921875, 'y9 l0+': 13567.1630859375, 'y13 l0+': 9119.26171875, 'y11 l0+': 15385.4248046875, 'y10 l0+': 7233.4287109375}, 23.329964: {'b4 h0+': 6373.9130859375, 'b8 h0+': 5655.77490234375, 'b7 h0+': 5740.88330078125, 'b5 h0+': 10147.6884765625, 'b3 h0+': 6264.6962890625, 'b6 h0+': 11031.962890625}, 23.240388: {'b4 h0+': 5602.08984375, 'b7 h0+': 8834.6865234375, 'b5 h0+': 11086.3310546875, 'b6 h0+': 8658.94921875, 'b2 h0+': 5011.15673828125}, 22.960066: {'y14 l0+': 10134.4189453125, 'y9 l0+': 34771.359375, 'y8 l0+': 7413.35791015625, 'y14 l1+': 3108.365478515625, 'y12 l0+': 26661.80078125, 'y9 l1+': 7776.69482421875, 'y16 l0+': 5923.12255859375, 'y7 l0+': 4948.7265625, 'y13 l1+': 9124.892578125, 'y11 l1+': 14329.333984375, 'y11 l0+': 44540.1796875, 'y13 l0+': 22150.560546875, 'y10 l0+': 20406.873046875, 'y12 l1+': 6211.6572265625}, 22.572267: {'y14 l0+': 36180.88671875, 'y9 l0+': 146540.171875, 'y5 l0+': 12797.978515625, 'y14 l1+': 15991.26953125, 'y9 l1+': 26751.947265625, 'y2 l0+': 7371.0712890625, 'y3 l0+': 16941.630859375, 'y13 l1+': 32177.794921875, 'y17 l1++': 13839.943359375, 'y11 l1+': 57276.640625, 'y11 l0+': 177767.28125, 'y4 l0+': 7410.94287109375, 'y10 l0+': 100166.1953125, 'y15 l1+': 4727.4140625, 'y12 l1+': 34872.67578125, 'y12 l0+': 91897.796875, 'y6 l0+': 5441.56201171875, 'y8 l0+': 20967.11328125, 'y10 l1+': 19234.068359375, 'y16 l0+': 18957.765625, 'y15 l0+': 16724.29296875, 'y7 l0+': 22175.62890625, 'y17 l0++': 8972.064453125, 'y13 l0+': 72237.5, 'y16 l1+': 7879.26513671875}, 23.371637: {'y9 l0+': 19239.091796875, 'y17 l2++': 32082.197265625, 'y11 l0+': 24289.78125, 'y12 l0+': 10206.919921875, 'y13 l0+': 10370.1533203125, 'y10 l0+': 14362.076171875, 'y17 l1++': 60461.6171875, 'y17 l1+': 9182.1396484375}, 22.634544: {'b5 h1+': 5663.4833984375, 'b7 h1+': 6386.38427734375, 'b6 h1+': 6486.1396484375, 'b5 h0+': 56455.984375, 'b9 h0+': 25551.291015625, 'b7 h0+': 40900.96875, 'b4 h0+': 28491.94140625, 'b8 h0+': 21744.76171875, 'b3 h0+': 6724.2216796875, 'b2 h0+': 23700.21875, 'b9 h1+': 5946.521484375, 'b8 h1+': 8667.509765625, 'b6 h0+': 39456.0234375}, 23.287187: {'b4 h0+': 7933.375, 'b2 h0+': 6749.25390625, 'b7 h0+': 9459.8505859375, 'b9 h0+': 6715.56982421875, 'b5 h0+': 14028.0078125, 'b6 h0+': 9750.0732421875}, 23.105453: {'y14 l0+': 5184.09912109375, 'y9 l0+': 22581.68359375, 'y8 l0+': 2804.0234375, 'y12 l0+': 13586.5693359375, 'y9 l1+': 4621.416015625, 'y16 l0+': 3031.02978515625, 'y7 l0+': 4784.71240234375, 'y13 l1+': 6763.78125, 'y11 l1+': 9236.150390625, 'y11 l0+': 32088.423828125, 'y13 l0+': 12844.7548828125, 'y10 l0+': 15707.26953125, 'y12 l1+': 6545.6337890625}, 22.519892: {'b5 h1+': 8510.1044921875, 'b7 h1+': 18370.86328125, 'b6 h1+': 10793.5341796875, 'b5 h0+': 118880.828125, 'b9 h0+': 59258.390625, 'b7 h0+': 104845.0625, 'b4 h0+': 55632.046875, 'b8 h0+': 61762.14453125, 'b3 h0+': 22855.349609375, 'y11 l1+': 7322.544921875, 'b2 h0+': 45708.26171875, 'b9 h1+': 15506.720703125, 'b8 h1+': 16203.345703125, 'b10 h0+': 10781.091796875, 'b6 h0+': 103506.1640625}, 23.993255: {'b4 h0+': 4893.1474609375, 'b9 h0+': 4833.73046875, 'b5 h0+': 8232.1875, 'b6 h0+': 8700.1533203125, 'b2 h0+': 5412.828125}}, 'MS1': {}}, 'SSEDPNEDIVER++': {'MS2': {11.379198: {'y9 l0+': 551251.75, 'y9 l0++': 225360.890625, 'y8 l1+': 347870.46875, 'y8 l0++': 1093021.25, 'y9 l1+': 112770.1640625, 'y10 l1+': 48113.69921875, 'y2 l0+': 156787.390625, 'y3 l0+': 244145.1875, 'y7 l0+': 180803.078125, 'y8 l0+': 2140608.5, 'y5 l0+': 254384.046875, 'y4 l0+': 100115.3671875, 'y6 l0+': 146867.296875, 'y10 l0+': 153993.21875, 'y8 l1++': 153571.3125}, 11.117339: {'y9 l0+': 311915.0, 'y5 l0+': 137951.359375, 'y8 l1+': 216298.59375, 'y5 l1+': 19669.083984375, 'y9 l1+': 76896.8125, 'y9 l0++': 114249.2421875, 'y2 l0+': 89602.5703125, 'y3 l1+': 7290.69482421875, 'y11 l0++': 16931.349609375, 'y3 l0+': 132501.796875, 'y6 l0+': 80285.890625, 'y4 l1+': 9733.4521484375, 'y10 l0+': 76138.6875, 'y8 l0+': 1062837.125, 'y8 l0++': 542676.25, 'y10 l1+': 20789.013671875, 'y9 l1++': 22910.54296875, 'y4 l0+': 77265.4140625, 'y10 l0++': 21695.13671875, 'y7 l1+': 16529.787109375, 'y7 l0+': 112596.7109375, 'y6 l1+': 13464.0322265625, 'y8 l1++': 102424.2578125}, 12.055869: {'y5 l2+': 4293.802734375, 'y6 l0+': 16104.90625, 'y6 l1+': 12115.2119140625, 'y6 l2+': 5125.267578125}, 11.211154: {'y10 l0++': 112025.734375, 'y9 l0+': 1526251.5, 'y8 l0+': 5473128.5, 'y8 l1+': 1068542.875, 'y5 l1+': 74474.7734375, 'y8 l0++': 2696013.25, 'y9 l1+': 323523.5, 'y9 l0++': 572964.125, 'y2 l0+': 488009.59375, 'y9 l1++': 127940.8203125, 'y3 l0+': 616830.4375, 'y7 l1+': 107917.4375, 'y7 l0+': 558316.875, 'y5 l0+': 662936.5625, 'y4 l0+': 396152.21875, 'y8 l1++': 560464.375, 'y6 l0+': 367222.25, 'y10 l0+': 350779.125, 'y10 l1+': 109015.84375}, 12.088127: {'y6 l0+': 6869.287109375}, 11.300961: {'b4 h0+': 1830040.25, 'b9 h0+': 167595.078125, 'b7 h0+': 229160.890625, 'b3 h0+': 610619.75, 'b2 h0+': 1859054.5}, 11.119583: {'b4 h0+': 551865.0625, 'b8 h0+': 24329.0390625, 'b4 h1+': 44379.51171875, 'b2 h0+': 541848.8125, 'b5 h0+': 40053.2109375, 'b10 h0+': 22903.119140625, 'b9 h0+': 48810.7109375, 'b7 h0+': 52218.265625, 'b3 h0+': 146122.203125}, 11.412908: {'y9 l0+': 218619.734375, 'y9 l0++': 79874.25, 'y8 l1+': 150948.03125, 'y8 l0++': 433738.15625, 'y9 l1+': 60482.359375, 'y10 l1+': 21362.267578125, 'y2 l0+': 74596.8203125, 'y9 l1++': 16253.6396484375, 'y3 l0+': 102713.6484375, 'y7 l0+': 85010.3671875, 'y8 l0+': 859961.6875, 'y5 l0+': 130139.8515625, 'y4 l0+': 40676.578125, 'y6 l0+': 61348.11328125, 'y10 l0+': 52391.05859375, 'y8 l1++': 71075.125}, 11.449535: {'y9 l0+': 73113.21875, 'y9 l0++': 24009.59375, 'y8 l1+': 41088.0078125, 'y8 l0++': 123902.71875, 'y9 l1++': 5729.8525390625, 'y9 l1+': 17870.3984375, 'y10 l1+': 5078.1298828125, 'y2 l0+': 14070.1904296875, 'y10 l0++': 4042.301025390625, 'y3 l0+': 29908.85546875, 'y7 l0+': 25400.953125, 'y8 l0+': 260240.515625, 'y5 l0+': 31988.31640625, 'y4 l0+': 16443.572265625, 'y6 l0+': 15670.3837890625, 'y10 l0+': 17220.779296875, 'y8 l1++': 17339.490234375}, 11.074094: {'b4 h0+': 166195.28125, 'b8 h0+': 6907.57421875, 'b4 h1+': 11513.2587890625, 'b2 h0+': 144453.03125, 'b7 h0+': 14693.427734375, 'b10 h0+': 5483.384765625, 'b9 h0+': 12589.818359375, 'b5 h0+': 6380.4892578125, 'b3 h0+': 38007.42578125, 'b6 h0+': 5866.19189453125}, 11.253128: {'y10 l0++': 139242.75, 'y9 l0+': 1804861.25, 'y8 l0+': 5923565.5, 'y9 l1+': 417348.0, 'y8 l1+': 1119235.125, 'y5 l1+': 62852.2265625, 'y8 l0++': 2900906.25, 'y10 l1+': 84322.15625, 'y9 l0++': 639199.0, 'y2 l0+': 452853.25, 'y9 l1++': 126859.9375, 'y3 l0+': 658765.5625, 'y7 l1+': 105424.5234375, 'y7 l0+': 655899.25, 'y11 l0++': 80611.6953125, 'y5 l0+': 686430.125, 'y4 l0+': 360982.5625, 'y6 l0+': 493947.03125, 'y10 l0+': 423678.6875, 'y8 l1++': 440412.4375}, 11.16295: {'y9 l0+': 874120.3125, 'y5 l0+': 367315.625, 'y8 l1+': 524231.75, 'y5 l1+': 57724.15625, 'y9 l1+': 207176.3125, 'y9 l0++': 285075.40625, 'y2 l0+': 227409.015625, 'y3 l0+': 345788.75, 'y11 l0++': 32679.951171875, 'y10 l1++': 29967.9140625, 'y4 l0+': 233300.3125, 'y6 l0+': 259057.078125, 'y10 l0+': 225258.96875, 'y8 l0+': 3171879.75, 'y8 l0++': 1556801.875, 'y8 l1++': 331553.9375, 'y9 l1++': 49397.59765625, 'y7 l0+': 384171.65625, 'y7 l1+': 60178.4765625, 'y10 l1+': 62443.5859375, 'y6 l1+': 39386.203125, 'y10 l0++': 40894.65625}, 11.341432: {'y10 l0++': 89828.375, 'y9 l0+': 1286164.125, 'y8 l0+': 4407595.5, 'y8 l1+': 852765.125, 'y5 l1+': 72621.296875, 'y8 l0++': 1889971.75, 'y9 l1+': 314136.71875, 'y9 l0++': 388042.75, 'y2 l0+': 333922.0, 'y9 l1++': 119503.28125, 'y3 l0+': 557429.25, 'y7 l1+': 70156.4921875, 'y7 l0+': 465385.90625, 'y5 l0+': 535408.5625, 'y4 l0+': 315410.0625, 'y8 l1++': 444176.0625, 'y6 l0+': 312849.6875, 'y10 l0+': 334220.84375, 'y10 l1+': 80225.2578125}, 11.485491: {'b4 h0+': 17331.3984375}, 12.125957: {'y6 l0+': 2647.498779296875}, 11.255374: {'b4 h0+': 1772873.25, 'b7 h0+': 188731.28125, 'b4 h1+': 142172.375, 'b9 h0+': 130897.7578125, 'b2 h0+': 1844218.375, 'b3 h0+': 454737.28125}, 11.451781: {'b4 h0+': 65468.67578125, 'b9 h0+': 6057.63916015625, 'b7 h0+': 6582.69482421875, 'b3 h0+': 14493.212890625, 'b2 h0+': 59696.7578125}, 11.298715: {'y9 l0+': 1495208.25, 'y8 l0+': 4956155.0, 'y8 l1+': 930318.375, 'y8 l0++': 2376688.0, 'y8 l1++': 435021.21875, 'y9 l0++': 578192.8125, 'y4 l0+': 313234.90625, 'y2 l0+': 406848.25, 'y9 l1++': 106463.4375, 'y7 l0+': 537282.875, 'y5 l0+': 701582.625, 'y6 l0+': 384725.53125, 'y3 l0+': 600752.5, 'y10 l0+': 369086.0625, 'y9 l1+': 241981.8125}, 12.238222: {'y5 l2+': 2461.18798828125}, 11.415154: {'b4 h0+': 221661.078125, 'b4 h1+': 15395.2685546875, 'b7 h0+': 17125.544921875, 'b10 h0+': 11580.140625, 'b9 h0+': 21418.255859375, 'b2 h0+': 208930.84375, 'b3 h0+': 59172.32421875}, 11.343677: {'b4 h0+': 1265458.0, 'b9 h0+': 103738.6328125, 'b7 h0+': 76802.3125, 'b3 h0+': 256043.796875, 'b2 h0+': 1140142.5}, 11.07185: {'y10 l0++': 7779.19189453125, 'y9 l0+': 114279.8671875, 'y8 l0+': 384373.28125, 'y8 l1+': 76235.625, 'y5 l1+': 9171.142578125, 'y8 l0++': 173647.6875, 'y9 l1+': 25120.939453125, 'y9 l0++': 36532.265625, 'y2 l0+': 35457.375, 'y9 l1++': 8283.068359375, 'y3 l0+': 40706.44140625, 'y7 l1+': 6975.08056640625, 'y7 l0+': 41705.859375, 'y5 l0+': 44677.16015625, 'y4 l0+': 27941.193359375, 'y8 l1++': 32993.16015625, 'y6 l0+': 30231.5390625, 'y10 l0+': 28970.884765625, 'y10 l1+': 7597.90673828125}, 11.483245: {'y9 l0+': 14483.6337890625, 'y7 l0+': 6021.5869140625, 'y5 l0+': 8805.697265625, 'y8 l0+': 61681.25, 'y8 l1+': 8524.6650390625, 'y6 l0+': 6499.59765625, 'y8 l0++': 27551.541015625, 'y3 l0+': 8133.1064453125, 'y9 l0++': 6951.86572265625, 'y2 l0+': 6699.99267578125}, 11.381444: {'b4 h0+': 615489.75, 'b7 h0+': 57737.52734375, 'b3 h0+': 185798.90625, 'b2 h0+': 680677.25}, 11.165195: {'b4 h0+': 1197730.125, 'b8 h0+': 61509.46875, 'b2 h0+': 1140154.0, 'b7 h0+': 103099.953125, 'b4 h1+': 95105.2734375, 'b9 h0+': 105129.515625, 'b5 h0+': 66890.4140625, 'b3 h0+': 287312.34375, 'b8 h1++': 49301.5, 'b6 h0+': 51303.02734375}, 11.213399: {'b4 h0+': 2016403.25, 'b8 h0+': 105794.109375, 'b7 h0+': 161939.375, 'b4 h1+': 170504.40625, 'b9 h0+': 167960.234375, 'b2 h0+': 2015227.25, 'b3 h0+': 555660.9375}}, 'MS1': {}}, 'SWPAVGNCSSALR++': {'MS2': {}, 'MS1': {}}, 'TVAACNLPIVR++': {'MS2': {11.268874: {'y9 l1++': 5103.326171875}}, 'MS1': {}}, 'FNAVLTNPQGDYDTSTGK++': {'MS2': {18.791707: {'y12 l2+': 8147.03173828125}, 18.318366: {'b4 h0+': 48093.84375, 'b2 h0+': 26928.529296875, 'b7 h0+': 19682.630859375, 'b4 h1+': 6417.17529296875, 'b5 h0+': 23917.0625, 'b3 h0+': 74380.234375, 'b6 h0+': 16948.958984375}, 18.210013: {'b5 h1+': 12718.3515625, 'b7 h1+': 27754.482421875, 'b2 h1+': 5959.50048828125, 'b5 h0+': 110587.140625, 'b4 h1+': 20573.3125, 'b9 h0+': 12486.498046875, 'b7 h0+': 110452.8046875, 'b16 h0+': 7423.009765625, 'b4 h0+': 211545.390625, 'b3 h1+': 19864.837890625, 'b6 h1+': 10887.5439453125, 'b11 h0+': 9223.1767578125, 'b2 h0+': 127055.4921875, 'b10 h0+': 6512.29638671875, 'b3 h0+': 344427.03125, 'b6 h0+': 73014.140625}, 18.147784: {'b5 h1+': 13464.572265625, 'b13 h0+': 9124.4658203125, 'b12 h0+': 5786.0166015625, 'b7 h1+': 16535.955078125, 'b6 h1+': 7304.94140625, 'b8 h0+': 9556.7587890625, 'b5 h0+': 80863.2734375, 'b4 h1+': 13812.90625, 'b9 h0+': 12665.7578125, 'b7 h0+': 84660.265625, 'b16 h0+': 7285.57568359375, 'b4 h0+': 146058.46875, 'b3 h1+': 19924.55859375, 'b11 h0+': 6520.26904296875, 'b2 h0+': 102005.3125, 'b3 h0+': 263760.09375, 'b6 h0+': 57563.52734375}, 18.372982: {'b3 h0+': 8359.158203125}, 18.757932: {'y12 l1++': 9515.8173828125}, 18.261473: {'b5 h1+': 10177.220703125, 'b7 h1+': 11958.177734375, 'b6 h1+': 5255.111328125, 'b8 h0+': 8778.5751953125, 'b5 h0+': 74141.3125, 'b4 h1+': 11947.8125, 'b9 h0+': 10170.708984375, 'b7 h0+': 72125.2421875, 'b14 h0+': 5203.4638671875, 'b4 h0+': 132758.421875, 'b3 h1+': 13435.275390625, 'b11 h0+': 6786.361328125, 'b2 h0+': 83702.015625, 'b10 h0+': 9650.89453125, 'b3 h0+': 229904.859375, 'b6 h0+': 46217.296875}, 18.082757: {'b4 h0+': 60102.4453125, 'b3 h1+': 8483.3251953125, 'b7 h1+': 5952.25537109375, 'b2 h0+': 41776.90234375, 'b7 h0+': 31046.83984375, 'b4 h1+': 6424.603515625, 'b5 h0+': 29249.93359375, 'b3 h0+': 106744.28125, 'b6 h0+': 26541.482421875}}, 'MS1': {}}}
    >>> compLib = create_compensationLibrary(peptideSettings, settings)
    A bunch of text...
    >>> CC = compLib['WSAGLTSSQVDLYIPK++'].compensationConstants
    >>> ref = {'y14 l2++': 1.2357433472284742, 'y10 l2+': 1.7023891343863502, 'y4 l1++': 1.7342556112536467, 'y12 h0++': 0.8323303392915962, 'y11 h1++': 1.1336821669869501, 'y2 h2+': 16.95414956508382, 'y6 h2+': 2.940889306061094, 'y10 l2++': 1.7016527571680717, 'y13 h2+': 1.3136703194966373, 'y5 h2+': 3.6622241853422413, 'y4 h1++': 1.7728802077998052, 'y12 l0++': 0.8259461766832151, 'y13 l0+': 0.8449120004342722, 'y7 l0++': 0.6525369769570856, 'y10 l0+': 0.7432740538589752, 'y12 l2++': 1.3667613774249254, 'y2 l2++': 16.491840694433076, 'y14 l0++': 0.8724448680658259, 'y13 h0+': 0.8498973200212949, 'y16 h1++': 1.0, 'y2 l2+': 16.523228223527376, 'y15 l0++': 0.9014749950082516, 'y13 l1+': 1.0679293654312407, 'y16 h0++': 1.0, 'y12 h2+': 1.3710324002348506, 'y13 l2++': 1.3081232818185902, 'y8 h0++': 0.7027475045128597, 'y5 h0+': 0.6007731683115595, 'y2 l0++': 0.46183467110790044, 'y6 h0++': 0.6300974429247793, 'y15 l2+': 1.1661727564764275, 'y11 h2++': 1.5465355298350725, 'y8 l0+': 0.692384054864364, 'y6 l1++': 1.405948650147822, 'y10 h0++': 0.7515827324167419, 'y8 l1++': 1.22567551804667, 'y5 h0++': 0.6009039381406013, 'y2 l0+': 0.46172979165633943, 'y15 h1++': 1.043034234012684, 'y7 l1++': 1.3051032859458325, 'y6 l1+': 1.4064080921943856, 'y4 h2++': 4.949637915192935, 'y10 l1++': 1.157585754577349, 'y11 l2+': 1.5391731034963723, 'y14 h0+': 0.8759515774578771, 'y4 h1+': 1.7738845250320894, 'y7 h0+': 0.664707404407357, 'y6 l2++': 2.9255216815063423, 'y3 h2++': 9.346178106337087, 'y9 l1++': 1.1889865888158906, 'y11 l2++': 1.5385688273206504, 'y4 l0+': 0.5487355369654731, 'y13 h2++': 1.3132050988811907, 'y3 h1+': 2.444840222393672, 'y4 h0+': 0.5625010450136524, 'y9 h2++': 1.8743429696179372, 'y15 l2++': 1.1658279898501545, 'y15 h2+': 1.1702857375323379, 'y5 h1++': 1.5371909659739473, 'y8 h2++': 2.070429923491537, 'y8 h1+': 1.2423081526035404, 'y12 l1++': 1.0802525438724173, 'y11 h0+': 0.7841217668040429, 'y16 h0+': 0.9998152787908181, 'y15 h0+': 0.9036525403772694, 'y8 l2++': 2.0557747803457436, 'y2 h2++': 16.921445097293827, 'y13 l1++': 1.0678158961882847, 'y8 h1++': 1.24201836012928, 'y10 h2++': 1.712723727191291, 'y14 h1++': 1.057432526667431, 'y16 h1+': 1.0000650154670363, 'y14 l2+': 1.236133810873526, 'y7 h1++': 1.322732690192651, 'y15 h2++': 1.1699191143222438, 'y3 h0++': 0.509015091863954, 'y3 l0+': 0.49540592284856744, 'y15 h0++': 0.9038254396264538, 'y13 l2+': 1.308562041381861, 'y11 l0+': 0.7769395228703242, 'y15 l1+': 1.038395217792486, 'y14 h1+': 1.0575466527289008, 'y4 l2+': 4.905611049527958, 'y2 h1+': 3.5320510722952236, 'y7 l1+': 1.3054478514383427, 'y10 h0+': 0.7514297388038317, 'y16 l2+': 1.0002494964831197, 'y9 l2+': 1.8625015615249094, 'y3 h0+': 0.5089004880452983, 'y16 l0++': 1.0, 'y5 l0+': 0.5872947347044313, 'y3 l2++': 9.134257428243412, 'y6 h0+': 0.6299622306127868, 'y7 l0+': 0.6523996320521745, 'y10 h2+': 1.7135022361584311, 'y3 l2+': 9.147513040852784, 'y11 l1++': 1.1226275963451835, 'y13 h1++': 1.0747358382368515, 'y14 l1++': 1.0516272033185345, 'y11 h1+': 1.1338639070148326, 'y2 h0+': 0.4753909581966258, 'y3 l0++': 0.4955172955268782, 'y4 l0++': 0.5488564183353715, 'y11 h2+': 1.5471743995296772, 'y15 l0+': 0.9013046146247351, 'y2 l1+': 3.3488163305191185, 'y14 l0+': 0.8722779508049894, 'y6 h2++': 2.9387993205714484, 'y13 l0++': 0.8450754898097119, 'y9 h2+': 1.8752702830972634, 'y16 l1++': 1.0, 'y8 l1+': 1.2259367254710718, 'y5 l1+': 1.5112372238579053, 'y4 h0++': 0.562625299165044, 'y4 l2++': 4.9007097799641235, 'y4 l1+': 1.735174040595471, 'y2 h1++': 3.52631153668765, 'y9 h1+': 1.204121463234435, 'y8 h0+': 0.7026012425093447, 'y12 h1+': 1.0876459488845187, 'y15 h1+': 1.0431353181977758, 'y9 h1++': 1.2038712751540799, 'y10 h1++': 1.1710397003304953, 'y9 l0+': 0.7175472382740405, 'y11 l0++': 0.777094112914068, 'y7 l2+': 2.4501431594851426, 'y5 l0++': 0.5874221093010794, 'y2 h0++': 0.4754990463376311, 'y11 h0++': 0.784279188389329, 'y7 l2++': 2.4486748614515665, 'y5 l2++': 3.6320070505118323, 'y16 h2+': 1.000265207348328, 'y6 h1++': 1.4285842052922928, 'y12 h2++': 1.3705251877819022, 'y16 l0+': 0.9998177965729644, 'y5 l2+': 3.634990357250837, 'y6 h1+': 1.4290867621331838, 'y8 h2+': 2.071552089207223, 'y16 h2++': 1.0, 'y4 h2+': 4.95476451617848, 'y3 h1++': 2.442460584127041, 'y13 h1+': 1.074865064308839, 'y12 h0+': 0.83216640202866, 'y3 l1++': 2.3540929231084653, 'y13 h0++': 0.850063524291167, 'y9 l2++': 1.8616234902055355, 'y6 l0++': 0.6172139998757586, 'y7 h2++': 2.4579187424098192, 'y16 l2++': 1.0, 'y14 h2++': 1.24068506950748, 'y5 l1++': 1.5106430241255997, 'y11 l1+': 1.1227889794841241, 'y2 l1++': 3.343824958232737, 'y12 l2+': 1.3672413756772155, 'y7 h0++': 0.6648480702293214, 'y14 l1+': 1.0517266825536096, 'y5 h2++': 3.659100422541519, 'y6 l0+': 0.6170821264171467, 'y16 l1+': 1.000054987820649, 'y8 l2+': 2.056838476714066, 'y9 l0++': 0.717693880291974, 'y7 h1+': 1.3231101865580965, 'y14 h2+': 1.241099970027311, 'y7 h2+': 2.4594540382923937, 'y5 h1+': 1.5378388381255013, 'y12 l0+': 0.8257850934127214, 'y10 h1+': 1.1712570444529298, 'y9 h0+': 0.72675925761864, 'y3 l1+': 2.356218759004883, 'y12 h1++': 1.0875052693206468, 'y9 h0++': 0.7269088854268524, 'y10 l0++': 0.7434241589995887, 'y6 l2+': 2.9275251251770733, 'y10 l1+': 1.1577796056617689, 'y14 h0++': 0.8761211109023641, 'y8 l0++': 0.692527233759618, 'y3 h2+': 9.360173229695084, 'y15 l1++': 1.0383077921191222, 'y9 l1+': 1.1892109241780022, 'y12 l1+': 1.0803770015327052}
    >>> all([str(CC[key])[:10]==str(ref[key])[:10] for key in CC])
    True
    >>> HI = compLib['WSAGLTSSQVDLYIPK++'].heavyIsoforms
    >>> ref = [{'S[33]': 0, 'O[16]': 0, 'S[32]': 0, 'C[13]': 0, 'O[17]': 0, 'H[1]': 0, 'C[12]': 0, 'N[14]': -2, 'S[34]': 0, 'O[18]': 0, 'N[15]': 2, 'H[2]': 0}, {'S[33]': 0, 'O[16]': 0, 'S[32]': -1, 'C[13]': 0, 'O[17]': 0, 'H[1]': 0, 'C[12]': 0, 'N[14]': -2, 'S[34]': 1, 'O[18]': 0, 'N[15]': 2, 'H[2]': 0}, {'S[33]': 1, 'O[16]': 0, 'S[32]': -1, 'C[13]': 0, 'O[17]': 0, 'H[1]': 0, 'C[12]': 0, 'N[14]': -2, 'S[34]': 0, 'O[18]': 0, 'N[15]': 2, 'H[2]': 0}, {'S[33]': 2, 'O[16]': 0, 'S[32]': -2, 'C[13]': 0, 'O[17]': 0, 'H[1]': 0, 'C[12]': 0, 'N[14]': -2, 'S[34]': 0, 'O[18]': 0, 'N[15]': 2, 'H[2]': 0}, {'S[33]': 0, 'O[16]': -1, 'S[32]': 0, 'C[13]': 0, 'O[17]': 0, 'H[1]': 0, 'C[12]': 0, 'N[14]': -2, 'S[34]': 0, 'O[18]': 1, 'N[15]': 2, 'H[2]': 0}, {'S[33]': 0, 'O[16]': -1, 'S[32]': 0, 'C[13]': 0, 'O[17]': 1, 'H[1]': 0, 'C[12]': 0, 'N[14]': -2, 'S[34]': 0, 'O[18]': 0, 'N[15]': 2, 'H[2]': 0}, {'S[33]': 1, 'O[16]': -1, 'S[32]': -1, 'C[13]': 0, 'O[17]': 1, 'H[1]': 0, 'C[12]': 0, 'N[14]': -2, 'S[34]': 0, 'O[18]': 0, 'N[15]': 2, 'H[2]': 0}, {'S[33]': 0, 'O[16]': -2, 'S[32]': 0, 'C[13]': 0, 'O[17]': 2, 'H[1]': 0, 'C[12]': 0, 'N[14]': -2, 'S[34]': 0, 'O[18]': 0, 'N[15]': 2, 'H[2]': 0}, {'S[33]': 0, 'O[16]': 0, 'S[32]': 0, 'C[13]': 0, 'O[17]': 0, 'H[1]': 0, 'C[12]': 0, 'N[14]': -3, 'S[34]': 0, 'O[18]': 0, 'N[15]': 3, 'H[2]': 0}, {'S[33]': 1, 'O[16]': 0, 'S[32]': -1, 'C[13]': 0, 'O[17]': 0, 'H[1]': 0, 'C[12]': 0, 'N[14]': -3, 'S[34]': 0, 'O[18]': 0, 'N[15]': 3, 'H[2]': 0}, {'S[33]': 0, 'O[16]': -1, 'S[32]': 0, 'C[13]': 0, 'O[17]': 1, 'H[1]': 0, 'C[12]': 0, 'N[14]': -3, 'S[34]': 0, 'O[18]': 0, 'N[15]': 3, 'H[2]': 0}, {'S[33]': 0, 'O[16]': 0, 'S[32]': 0, 'C[13]': 0, 'O[17]': 0, 'H[1]': 0, 'C[12]': 0, 'N[14]': -4, 'S[34]': 0, 'O[18]': 0, 'N[15]': 4, 'H[2]': 0}, {'S[33]': 0, 'O[16]': 0, 'S[32]': 0, 'C[13]': 0, 'O[17]': 0, 'H[1]': -1, 'C[12]': 0, 'N[14]': -2, 'S[34]': 0, 'O[18]': 0, 'N[15]': 2, 'H[2]': 1}, {'S[33]': 1, 'O[16]': 0, 'S[32]': -1, 'C[13]': 0, 'O[17]': 0, 'H[1]': -1, 'C[12]': 0, 'N[14]': -2, 'S[34]': 0, 'O[18]': 0, 'N[15]': 2, 'H[2]': 1}, {'S[33]': 0, 'O[16]': -1, 'S[32]': 0, 'C[13]': 0, 'O[17]': 1, 'H[1]': -1, 'C[12]': 0, 'N[14]': -2, 'S[34]': 0, 'O[18]': 0, 'N[15]': 2, 'H[2]': 1}, {'S[33]': 0, 'O[16]': 0, 'S[32]': 0, 'C[13]': 0, 'O[17]': 0, 'H[1]': -1, 'C[12]': 0, 'N[14]': -3, 'S[34]': 0, 'O[18]': 0, 'N[15]': 3, 'H[2]': 1}, {'S[33]': 0, 'O[16]': 0, 'S[32]': 0, 'C[13]': 0, 'O[17]': 0, 'H[1]': -2, 'C[12]': 0, 'N[14]': -2, 'S[34]': 0, 'O[18]': 0, 'N[15]': 2, 'H[2]': 2}, {'S[33]': 0, 'O[16]': 0, 'S[32]': 0, 'C[13]': 1, 'O[17]': 0, 'H[1]': 0, 'C[12]': -1, 'N[14]': -2, 'S[34]': 0, 'O[18]': 0, 'N[15]': 2, 'H[2]': 0}, {'S[33]': 1, 'O[16]': 0, 'S[32]': -1, 'C[13]': 1, 'O[17]': 0, 'H[1]': 0, 'C[12]': -1, 'N[14]': -2, 'S[34]': 0, 'O[18]': 0, 'N[15]': 2, 'H[2]': 0}, {'S[33]': 0, 'O[16]': -1, 'S[32]': 0, 'C[13]': 1, 'O[17]': 1, 'H[1]': 0, 'C[12]': -1, 'N[14]': -2, 'S[34]': 0, 'O[18]': 0, 'N[15]': 2, 'H[2]': 0}, {'S[33]': 0, 'O[16]': 0, 'S[32]': 0, 'C[13]': 1, 'O[17]': 0, 'H[1]': 0, 'C[12]': -1, 'N[14]': -3, 'S[34]': 0, 'O[18]': 0, 'N[15]': 3, 'H[2]': 0}, {'S[33]': 0, 'O[16]': 0, 'S[32]': 0, 'C[13]': 1, 'O[17]': 0, 'H[1]': -1, 'C[12]': -1, 'N[14]': -2, 'S[34]': 0, 'O[18]': 0, 'N[15]': 2, 'H[2]': 1}, {'S[33]': 0, 'O[16]': 0, 'S[32]': 0, 'C[13]': 2, 'O[17]': 0, 'H[1]': 0, 'C[12]': -2, 'N[14]': -2, 'S[34]': 0, 'O[18]': 0, 'N[15]': 2, 'H[2]': 0}]
    >>> all([str(HI[i][key])[:10]==str(ref[i][key])[:10] for i in range(len(HI)) for key in HI[i]])
    True
    >>> str(compLib['WSAGLTSSQVDLYIPK++'].heavyPrecursorDistribution[0])[:10]
    '0.41804198'
    >>> str(compLib['WSAGLTSSQVDLYIPK++'].heavyPrecursorDistribution[1])[:10]
    '0.38507346'
    >>> str(compLib['WSAGLTSSQVDLYIPK++'].heavyPrecursorDistribution[2])[:10]
    '0.19688454'
    >>> compLib['WSAGLTSSQVDLYIPK++'].includedPeaks
    ['l0', 'h0', 'l1', 'h1', 'l2', 'h2']
    >>> LI = compLib['WSAGLTSSQVDLYIPK++'].lightIsoforms
    >>> ref = [{'S[33]': 0, 'O[16]': 0, 'S[32]': 0, 'C[13]': 0, 'O[17]': 0, 'H[1]': 0, 'C[12]': 0, 'N[14]': 0, 'S[34]': 0, 'O[18]': 0, 'N[15]': 0, 'H[2]': 0}, {'S[33]': 0, 'O[16]': 0, 'S[32]': -1, 'C[13]': 0, 'O[17]': 0, 'H[1]': 0, 'C[12]': 0, 'N[14]': 0, 'S[34]': 1, 'O[18]': 0, 'N[15]': 0, 'H[2]': 0}, {'S[33]': 1, 'O[16]': 0, 'S[32]': -1, 'C[13]': 0, 'O[17]': 0, 'H[1]': 0, 'C[12]': 0, 'N[14]': 0, 'S[34]': 0, 'O[18]': 0, 'N[15]': 0, 'H[2]': 0}, {'S[33]': 2, 'O[16]': 0, 'S[32]': -2, 'C[13]': 0, 'O[17]': 0, 'H[1]': 0, 'C[12]': 0, 'N[14]': 0, 'S[34]': 0, 'O[18]': 0, 'N[15]': 0, 'H[2]': 0}, {'S[33]': 0, 'O[16]': -1, 'S[32]': 0, 'C[13]': 0, 'O[17]': 0, 'H[1]': 0, 'C[12]': 0, 'N[14]': 0, 'S[34]': 0, 'O[18]': 1, 'N[15]': 0, 'H[2]': 0}, {'S[33]': 0, 'O[16]': -1, 'S[32]': 0, 'C[13]': 0, 'O[17]': 1, 'H[1]': 0, 'C[12]': 0, 'N[14]': 0, 'S[34]': 0, 'O[18]': 0, 'N[15]': 0, 'H[2]': 0}, {'S[33]': 1, 'O[16]': -1, 'S[32]': -1, 'C[13]': 0, 'O[17]': 1, 'H[1]': 0, 'C[12]': 0, 'N[14]': 0, 'S[34]': 0, 'O[18]': 0, 'N[15]': 0, 'H[2]': 0}, {'S[33]': 0, 'O[16]': -2, 'S[32]': 0, 'C[13]': 0, 'O[17]': 2, 'H[1]': 0, 'C[12]': 0, 'N[14]': 0, 'S[34]': 0, 'O[18]': 0, 'N[15]': 0, 'H[2]': 0}, {'S[33]': 0, 'O[16]': 0, 'S[32]': 0, 'C[13]': 0, 'O[17]': 0, 'H[1]': 0, 'C[12]': 0, 'N[14]': -1, 'S[34]': 0, 'O[18]': 0, 'N[15]': 1, 'H[2]': 0}, {'S[33]': 1, 'O[16]': 0, 'S[32]': -1, 'C[13]': 0, 'O[17]': 0, 'H[1]': 0, 'C[12]': 0, 'N[14]': -1, 'S[34]': 0, 'O[18]': 0, 'N[15]': 1, 'H[2]': 0}, {'S[33]': 0, 'O[16]': -1, 'S[32]': 0, 'C[13]': 0, 'O[17]': 1, 'H[1]': 0, 'C[12]': 0, 'N[14]': -1, 'S[34]': 0, 'O[18]': 0, 'N[15]': 1, 'H[2]': 0}, {'S[33]': 0, 'O[16]': 0, 'S[32]': 0, 'C[13]': 0, 'O[17]': 0, 'H[1]': 0, 'C[12]': 0, 'N[14]': -2, 'S[34]': 0, 'O[18]': 0, 'N[15]': 2, 'H[2]': 0}, {'S[33]': 0, 'O[16]': 0, 'S[32]': 0, 'C[13]': 0, 'O[17]': 0, 'H[1]': -1, 'C[12]': 0, 'N[14]': 0, 'S[34]': 0, 'O[18]': 0, 'N[15]': 0, 'H[2]': 1}, {'S[33]': 1, 'O[16]': 0, 'S[32]': -1, 'C[13]': 0, 'O[17]': 0, 'H[1]': -1, 'C[12]': 0, 'N[14]': 0, 'S[34]': 0, 'O[18]': 0, 'N[15]': 0, 'H[2]': 1}, {'S[33]': 0, 'O[16]': -1, 'S[32]': 0, 'C[13]': 0, 'O[17]': 1, 'H[1]': -1, 'C[12]': 0, 'N[14]': 0, 'S[34]': 0, 'O[18]': 0, 'N[15]': 0, 'H[2]': 1}, {'S[33]': 0, 'O[16]': 0, 'S[32]': 0, 'C[13]': 0, 'O[17]': 0, 'H[1]': -1, 'C[12]': 0, 'N[14]': -1, 'S[34]': 0, 'O[18]': 0, 'N[15]': 1, 'H[2]': 1}, {'S[33]': 0, 'O[16]': 0, 'S[32]': 0, 'C[13]': 0, 'O[17]': 0, 'H[1]': -2, 'C[12]': 0, 'N[14]': 0, 'S[34]': 0, 'O[18]': 0, 'N[15]': 0, 'H[2]': 2}, {'S[33]': 0, 'O[16]': 0, 'S[32]': 0, 'C[13]': 1, 'O[17]': 0, 'H[1]': 0, 'C[12]': -1, 'N[14]': 0, 'S[34]': 0, 'O[18]': 0, 'N[15]': 0, 'H[2]': 0}, {'S[33]': 1, 'O[16]': 0, 'S[32]': -1, 'C[13]': 1, 'O[17]': 0, 'H[1]': 0, 'C[12]': -1, 'N[14]': 0, 'S[34]': 0, 'O[18]': 0, 'N[15]': 0, 'H[2]': 0}, {'S[33]': 0, 'O[16]': -1, 'S[32]': 0, 'C[13]': 1, 'O[17]': 1, 'H[1]': 0, 'C[12]': -1, 'N[14]': 0, 'S[34]': 0, 'O[18]': 0, 'N[15]': 0, 'H[2]': 0}, {'S[33]': 0, 'O[16]': 0, 'S[32]': 0, 'C[13]': 1, 'O[17]': 0, 'H[1]': 0, 'C[12]': -1, 'N[14]': -1, 'S[34]': 0, 'O[18]': 0, 'N[15]': 1, 'H[2]': 0}, {'S[33]': 0, 'O[16]': 0, 'S[32]': 0, 'C[13]': 1, 'O[17]': 0, 'H[1]': -1, 'C[12]': -1, 'N[14]': 0, 'S[34]': 0, 'O[18]': 0, 'N[15]': 0, 'H[2]': 1}, {'S[33]': 0, 'O[16]': 0, 'S[32]': 0, 'C[13]': 2, 'O[17]': 0, 'H[1]': 0, 'C[12]': -2, 'N[14]': 0, 'S[34]': 0, 'O[18]': 0, 'N[15]': 0, 'H[2]': 0}]
    >>> all([str(LI[i][key])[:10]==str(ref[i][key])[:10] for i in range(len(LI)) for key in LI[i]])
    True
    >>> compLib['WSAGLTSSQVDLYIPK++'].overlappingPeaks=={'y14 l2++': 'y14 h0++', 'y10 l2+': 'y10 h0+', 'y5 l0++': None, 'y12 h0++': 'y12 l2++', 'y11 h1++': None, 'y2 h2+': None, 'y6 h2+': None, 'y10 l2++': 'y10 h0++', 'y13 h2+': None, 'y5 h2+': None, 'y7 l1++': None, 'y7 h0++': 'y7 l2++', 'y12 l0++': None, 'y13 l0+': None, 'y7 l0++': None, 'y10 l0+': None, 'y12 l2++': 'y12 h0++', 'y3 h0++': 'y3 l2++', 'y7 l0+': None, 'y13 h0+': 'y13 l2+', 'y16 h1++': None, 'y2 l2+': 'y2 h0+', 'y11 h1+': None, 'y13 l1+': None, 'y16 h0++': 'y16 l2++', 'y12 h2+': None, 'y13 l2++': 'y13 h0++', 'y8 h0++': 'y8 l2++', 'y5 h0+': 'y5 l2+', 'y2 l0++': None, 'y6 h0++': 'y6 l2++', 'y15 l2+': 'y15 h0+', 'y8 l0+': None, 'l2': 'h0', 'y7 h1+': None, 'y10 h0++': 'y10 l2++', 'y8 l1++': None, 'y5 h0++': 'y5 l2++', 'y2 l0+': None, 'y14 h2+': None, 'y15 h1++': None, 'y4 h1++': None, 'y6 l1+': None, 'y4 h2++': None, 'y10 l1++': None, 'y3 h1+': None, 'y14 h0+': 'y14 l2+', 'y4 h1+': None, 'y7 h0+': 'y7 l2+', 'y6 l2++': 'y6 h0++', 'y3 h2++': None, 'y9 l1++': None, 'h0': 'l2', 'y11 l2++': 'y11 h0++', 'y4 l0+': None, 'y15 h2++': None, 'y11 l2+': 'y11 h0+', 'y4 h0+': 'y4 l2+', 'y9 h2++': None, 'y15 l2++': 'y15 h0++', 'y15 h2+': None, 'y5 h1++': None, 'y9 h1+': None, 'y12 l1++': None, 'y11 h0+': 'y11 l2+', 'y16 h0+': 'y16 l2+', 'y15 h0+': 'y15 l2+', 'y8 l2++': 'y8 h0++', 'y2 h2++': None, 'y13 l1++': None, 'y8 h1++': None, 'y10 h2++': None, 'y14 h1++': None, 'y16 h1+': None, 'y14 l2+': 'y14 h0+', 'y7 h1++': None, 'y13 h2++': None, 'y2 l2++': 'y2 h0++', 'y3 l0+': None, 'y15 h0++': 'y15 l2++', 'y7 l1+': None, 'y14 l0++': None, 'y11 l0+': None, 'y15 l1+': None, 'y14 h1+': None, 'y4 l2+': 'y4 h0+', 'y2 h1+': None, 'y13 l2+': 'y13 h0+', 'y10 h0+': 'y10 l2+', 'y16 l2+': 'y16 h0+', 'y9 l2+': 'y9 h0+', 'y14 l0+': None, 'l0': None, 'y16 l0++': None, 'y5 l0+': None, 'y3 l2++': 'y3 h0++', 'y6 h0+': 'y6 l2+', 'y10 h2+': None, 'y3 l2+': 'y3 h0+', 'y11 l1++': None, 'y13 h1++': None, 'y14 l1++': None, 'y15 l0++': None, 'y2 h0+': 'y2 l2+', 'y3 l0++': None, 'y4 l0++': None, 'y11 h2+': None, 'y15 l0+': None, 'y2 l1+': None, 'y8 h1+': None, 'y6 h2++': None, 'y13 l0++': None, 'y9 h2+': None, 'y16 l1++': None, 'y8 l1+': None, 'y5 l1+': None, 'y4 h0++': 'y4 l2++', 'y5 l1++': None, 'y4 l1+': None, 'y2 h1++': None, 'y8 h0+': 'y8 l2+', 'y12 h1+': None, 'y9 h1++': None, 'y5 h2++': None, 'y9 l0+': None, 'y11 l0++': None, 'y7 l2+': 'y7 h0+', 'y4 l1++': None, 'y2 h0++': 'y2 l2++', 'y11 h0++': 'y11 l2++', 'y7 l2++': 'y7 h0++', 'y5 l2++': 'y5 h0++', 'y16 h2+': None, 'y8 h2++': None, 'y12 h2++': None, 'y16 l0+': None, 'y5 l2+': 'y5 h0+', 'y6 h1+': None, 'y9 h0++': 'y9 l2++', 'y16 h2++': None, 'y4 h2+': None, 'h1': None, 'y3 h1++': None, 'y13 h1+': None, 'y12 h0+': 'y12 l2+', 'y3 l1++': None, 'y13 h0++': 'y13 l2++', 'y9 l2++': 'y9 h0++', 'y6 l0++': None, 'y7 h2++': None, 'y16 l2++': 'y16 h0++', 'y14 h2++': None, 'y11 l1+': None, 'y2 l1++': None, 'y12 l2+': 'y12 h0+', 'y15 h1+': None, 'y14 l1+': None, 'y6 h1++': None, 'y10 h1++': None, 'y6 l0+': None, 'y16 l1+': None, 'y8 l2+': 'y8 h0+', 'y9 l0++': None, 'y6 l1++': None, 'y3 h0+': 'y3 l2+', 'y7 h2+': None, 'y8 h2+': None, 'y5 h1+': None, 'y12 l0+': None, 'y10 h1+': None, 'y9 h0+': 'y9 l2+', 'y3 l1+': None, 'y12 h1++': None, 'l1': None, 'y10 l0++': None, 'y11 h2++': None, 'y10 l1+': None, 'y14 h0++': 'y14 l2++', 'y8 l0++': None, 'y3 h2+': None, 'y15 l1++': None, 'y9 l1+': None, 'y12 l1+': None, 'y4 l2++': 'y4 h0++', 'h2': None, 'y6 l2+': 'y6 h0+'}
    True
    >>> compLib['WSAGLTSSQVDLYIPK++'].precursorDistribution
    [0.40214660850390893, 0.39002048124937416, 0.20783291024671688]
    >>> CC = compLib['FVYHLSDLCK+++'].compensationConstants
    >>> ref = {'y10 l2+': 1.0006128338413063, 'y4 l1++': 1.794607739993808, 'y9 l0+': 0.9155393758452934, 'y2 h0+': 0.6030512457732322, 'y2 h2+': 2.9922645823280063, 'y6 h2+': 1.7841296511665719, 'y3 l0++': 0.6331228418240626, 'y4 l0++': 0.6655519293428189, 'y5 h2+': 2.047802031504586, 'y4 h1++': 1.8472407675843339, 'y6 l0++': 0.737772358439925, 'y2 l1+': 2.8980465679765897, 'y8 h1+': 1.134681460687604, 'y6 h2++': 1.783372409366896, 'y10 l0+': 0.9996140942429314, 'y3 h0++': 0.6439781291804547, 'y7 l0+': 0.7916385142833069, 'y4 l2+': 2.3071597669040376, 'y2 l2+': 3.0950655411361825, 'y8 l1+': 1.1238951574220142, 'y5 l2++': 2.0763097761905462, 'y4 h0++': 0.6755291084855758, 'y5 l1++': 1.632437228453558, 'y4 l1+': 1.7957836318364997, 'y2 l1++': 2.894111476323344, 'y5 h0+': 0.6997201430888307, 'y2 h1++': 3.05401152776191, 'y2 l0++': 0.5917133634169215, 'y6 h0++': 0.7457597802730835, 'y8 h0+': 0.8715057917793761, 'y3 l1++': 2.091209600863749, 'y8 l0+': 0.868890989498753, 'y9 h1++': 1.0760580361171663, 'y5 h2++': 2.0469555972622597, 'y8 h0++': 0.8716823186345406, 'y7 h1+': 1.2827340273662993, 'y10 h0++': 0.9998052547427604, 'y8 l1++': 1.1236529993942648, 'y5 h0++': 0.6998690554054698, 'y10 h0+': 0.9996105361919443, 'y7 l2+': 1.549264723703335, 'y5 l0++': 0.6909940834023125, 'y2 l0+': 0.5915843568688826, 'y7 l2++': 1.5486319465067735, 'y7 l1++': 1.2578701799858458, 'y6 l1+': 1.4150458025722732, 'y4 h2++': 2.2640404254350917, 'y6 h1++': 1.4414657927404075, 'y10 l1++': 1.0001292075131618, 'y3 h1+': 2.1665565460227345, 'y4 h1+': 1.8485608238932, 'y7 h0+': 0.7952233368469389, 'y6 l2++': 1.7997269932518873, 'y3 h2++': 2.6320887548946863, 'y8 h2+': 1.2831451872513224, 'y9 h2+': 1.1641183196515243, 'y8 h2++': 1.2826606743411042, 'y4 h2+': 2.264978532013838, 'y3 h1++': 2.1645530453609427, 'y4 l0+': 0.6654095441495428, 'y7 l0++': 0.7918016913583837, 'y9 h0++': 0.916914112218578, 'y4 h0+': 0.6753842599394696, 'y9 h2++': 1.1637048680375945, 'y9 l1+': 1.069632412428205, 'y9 l2++': 1.1611880833695547, 'y5 h1++': 1.676910208795631, 'y7 h2++': 1.5525254702624058, 'y9 h1+': 1.076271958360603, 'y8 l0++': 0.8690657397529351, 'y5 l2+': 2.0771875704513656, 'y8 h1++': 1.134409859685743, 'y2 h2++': 2.9914152171138375, 'y3 l2++': 2.69372558124503, 'y6 h1+': 1.4421090913644956, 'y2 h0++': 0.6031828637685044, 'y7 h0++': 0.7953882121829775, 'y10 h2++': 1.0003184509883332, 'y10 h1++': 1.0001448049051973, 'y9 l1++': 1.069441871234485, 'y6 l0+': 0.7376178552160957, 'y7 h1++': 1.2822977524612873, 'y2 l2++': 3.0940900403239207, 'y8 l2+': 1.2810299387516872, 'y3 l0+': 0.6329859322669239, 'y9 l0++': 0.9157207854834086, 'y6 l1++': 1.4144680922346067, 'y8 l2++': 1.2805587762160378, 'y7 l1+': 1.258254858559398, 'y3 h0+': 0.6438386443535574, 'y7 h2+': 1.5531682454250062, 'y3 l2+': 2.6948748013965185, 'y5 h1+': 1.6779183225448422, 'y2 h1+': 3.0586027037213177, 'y10 h1+': 1.0002897349245914, 'y9 h0+': 0.916731000960836, 'y3 l1+': 2.0929843302788025, 'y5 l1+': 1.6333338917066968, 'y9 l2+': 1.1615875599572645, 'y10 l0++': 0.9998070341722968, 'y6 l2+': 1.800493330400712, 'y10 l1+': 1.0002585242386994, 'y5 l0+': 0.690847483187583, 'y3 h2+': 2.63316034690671, 'y10 l2++': 1.0003063314463652, 'y6 h0+': 0.7456030011248194, 'y10 h2+': 1.0006370812881984, 'y4 l2++': 2.306173866562879}
    >>> all([str(CC[key])[:10]==str(ref[key])[:10] for key in CC])
    True
    >>> HI = compLib['FVYHLSDLCK+++'].heavyIsoforms
    >>> ref = [{'S[33]': 0, 'O[16]': 0, 'S[32]': 0, 'C[13]': 0, 'O[17]': 0, 'H[1]': 0, 'C[12]': 0, 'N[14]': -2, 'S[34]': 0, 'O[18]': 0, 'N[15]': 2, 'H[2]': 0}, {'S[33]': 0, 'O[16]': 0, 'S[32]': -1, 'C[13]': 0, 'O[17]': 0, 'H[1]': 0, 'C[12]': 0, 'N[14]': -2, 'S[34]': 1, 'O[18]': 0, 'N[15]': 2, 'H[2]': 0}, {'S[33]': 1, 'O[16]': 0, 'S[32]': -1, 'C[13]': 0, 'O[17]': 0, 'H[1]': 0, 'C[12]': 0, 'N[14]': -2, 'S[34]': 0, 'O[18]': 0, 'N[15]': 2, 'H[2]': 0}, {'S[33]': 2, 'O[16]': 0, 'S[32]': -2, 'C[13]': 0, 'O[17]': 0, 'H[1]': 0, 'C[12]': 0, 'N[14]': -2, 'S[34]': 0, 'O[18]': 0, 'N[15]': 2, 'H[2]': 0}, {'S[33]': 0, 'O[16]': -1, 'S[32]': 0, 'C[13]': 0, 'O[17]': 0, 'H[1]': 0, 'C[12]': 0, 'N[14]': -2, 'S[34]': 0, 'O[18]': 1, 'N[15]': 2, 'H[2]': 0}, {'S[33]': 0, 'O[16]': -1, 'S[32]': 0, 'C[13]': 0, 'O[17]': 1, 'H[1]': 0, 'C[12]': 0, 'N[14]': -2, 'S[34]': 0, 'O[18]': 0, 'N[15]': 2, 'H[2]': 0}, {'S[33]': 1, 'O[16]': -1, 'S[32]': -1, 'C[13]': 0, 'O[17]': 1, 'H[1]': 0, 'C[12]': 0, 'N[14]': -2, 'S[34]': 0, 'O[18]': 0, 'N[15]': 2, 'H[2]': 0}, {'S[33]': 0, 'O[16]': -2, 'S[32]': 0, 'C[13]': 0, 'O[17]': 2, 'H[1]': 0, 'C[12]': 0, 'N[14]': -2, 'S[34]': 0, 'O[18]': 0, 'N[15]': 2, 'H[2]': 0}, {'S[33]': 0, 'O[16]': 0, 'S[32]': 0, 'C[13]': 0, 'O[17]': 0, 'H[1]': 0, 'C[12]': 0, 'N[14]': -3, 'S[34]': 0, 'O[18]': 0, 'N[15]': 3, 'H[2]': 0}, {'S[33]': 1, 'O[16]': 0, 'S[32]': -1, 'C[13]': 0, 'O[17]': 0, 'H[1]': 0, 'C[12]': 0, 'N[14]': -3, 'S[34]': 0, 'O[18]': 0, 'N[15]': 3, 'H[2]': 0}, {'S[33]': 0, 'O[16]': -1, 'S[32]': 0, 'C[13]': 0, 'O[17]': 1, 'H[1]': 0, 'C[12]': 0, 'N[14]': -3, 'S[34]': 0, 'O[18]': 0, 'N[15]': 3, 'H[2]': 0}, {'S[33]': 0, 'O[16]': 0, 'S[32]': 0, 'C[13]': 0, 'O[17]': 0, 'H[1]': 0, 'C[12]': 0, 'N[14]': -4, 'S[34]': 0, 'O[18]': 0, 'N[15]': 4, 'H[2]': 0}, {'S[33]': 0, 'O[16]': 0, 'S[32]': 0, 'C[13]': 0, 'O[17]': 0, 'H[1]': -1, 'C[12]': 0, 'N[14]': -2, 'S[34]': 0, 'O[18]': 0, 'N[15]': 2, 'H[2]': 1}, {'S[33]': 1, 'O[16]': 0, 'S[32]': -1, 'C[13]': 0, 'O[17]': 0, 'H[1]': -1, 'C[12]': 0, 'N[14]': -2, 'S[34]': 0, 'O[18]': 0, 'N[15]': 2, 'H[2]': 1}, {'S[33]': 0, 'O[16]': -1, 'S[32]': 0, 'C[13]': 0, 'O[17]': 1, 'H[1]': -1, 'C[12]': 0, 'N[14]': -2, 'S[34]': 0, 'O[18]': 0, 'N[15]': 2, 'H[2]': 1}, {'S[33]': 0, 'O[16]': 0, 'S[32]': 0, 'C[13]': 0, 'O[17]': 0, 'H[1]': -1, 'C[12]': 0, 'N[14]': -3, 'S[34]': 0, 'O[18]': 0, 'N[15]': 3, 'H[2]': 1}, {'S[33]': 0, 'O[16]': 0, 'S[32]': 0, 'C[13]': 0, 'O[17]': 0, 'H[1]': -2, 'C[12]': 0, 'N[14]': -2, 'S[34]': 0, 'O[18]': 0, 'N[15]': 2, 'H[2]': 2}, {'S[33]': 0, 'O[16]': 0, 'S[32]': 0, 'C[13]': 1, 'O[17]': 0, 'H[1]': 0, 'C[12]': -1, 'N[14]': -2, 'S[34]': 0, 'O[18]': 0, 'N[15]': 2, 'H[2]': 0}, {'S[33]': 1, 'O[16]': 0, 'S[32]': -1, 'C[13]': 1, 'O[17]': 0, 'H[1]': 0, 'C[12]': -1, 'N[14]': -2, 'S[34]': 0, 'O[18]': 0, 'N[15]': 2, 'H[2]': 0}, {'S[33]': 0, 'O[16]': -1, 'S[32]': 0, 'C[13]': 1, 'O[17]': 1, 'H[1]': 0, 'C[12]': -1, 'N[14]': -2, 'S[34]': 0, 'O[18]': 0, 'N[15]': 2, 'H[2]': 0}, {'S[33]': 0, 'O[16]': 0, 'S[32]': 0, 'C[13]': 1, 'O[17]': 0, 'H[1]': 0, 'C[12]': -1, 'N[14]': -3, 'S[34]': 0, 'O[18]': 0, 'N[15]': 3, 'H[2]': 0}, {'S[33]': 0, 'O[16]': 0, 'S[32]': 0, 'C[13]': 1, 'O[17]': 0, 'H[1]': -1, 'C[12]': -1, 'N[14]': -2, 'S[34]': 0, 'O[18]': 0, 'N[15]': 2, 'H[2]': 1}, {'S[33]': 0, 'O[16]': 0, 'S[32]': 0, 'C[13]': 2, 'O[17]': 0, 'H[1]': 0, 'C[12]': -2, 'N[14]': -2, 'S[34]': 0, 'O[18]': 0, 'N[15]': 2, 'H[2]': 0}]
    >>> all([str(HI[i][key])[:10]==str(ref[i][key])[:10] for i in range(len(HI)) for key in HI[i]])
    True
    >>> compLib['FVYHLSDLCK+++'].heavyPrecursorDistribution
    [0.5047312187058874, 0.34197699455207126, 0.15329178674204127]
    >>> compLib['FVYHLSDLCK+++'].includedPeaks
    ['l0', 'h0', 'l1', 'h1', 'l2', 'h2']
    >>> LI = compLib['FVYHLSDLCK+++'].lightIsoforms
    >>> ref = [{'S[33]': 0, 'O[16]': 0, 'S[32]': 0, 'C[13]': 0, 'O[17]': 0, 'H[1]': 0, 'C[12]': 0, 'N[14]': 0, 'S[34]': 0, 'O[18]': 0, 'N[15]': 0, 'H[2]': 0}, {'S[33]': 0, 'O[16]': 0, 'S[32]': -1, 'C[13]': 0, 'O[17]': 0, 'H[1]': 0, 'C[12]': 0, 'N[14]': 0, 'S[34]': 1, 'O[18]': 0, 'N[15]': 0, 'H[2]': 0}, {'S[33]': 1, 'O[16]': 0, 'S[32]': -1, 'C[13]': 0, 'O[17]': 0, 'H[1]': 0, 'C[12]': 0, 'N[14]': 0, 'S[34]': 0, 'O[18]': 0, 'N[15]': 0, 'H[2]': 0}, {'S[33]': 2, 'O[16]': 0, 'S[32]': -2, 'C[13]': 0, 'O[17]': 0, 'H[1]': 0, 'C[12]': 0, 'N[14]': 0, 'S[34]': 0, 'O[18]': 0, 'N[15]': 0, 'H[2]': 0}, {'S[33]': 0, 'O[16]': -1, 'S[32]': 0, 'C[13]': 0, 'O[17]': 0, 'H[1]': 0, 'C[12]': 0, 'N[14]': 0, 'S[34]': 0, 'O[18]': 1, 'N[15]': 0, 'H[2]': 0}, {'S[33]': 0, 'O[16]': -1, 'S[32]': 0, 'C[13]': 0, 'O[17]': 1, 'H[1]': 0, 'C[12]': 0, 'N[14]': 0, 'S[34]': 0, 'O[18]': 0, 'N[15]': 0, 'H[2]': 0}, {'S[33]': 1, 'O[16]': -1, 'S[32]': -1, 'C[13]': 0, 'O[17]': 1, 'H[1]': 0, 'C[12]': 0, 'N[14]': 0, 'S[34]': 0, 'O[18]': 0, 'N[15]': 0, 'H[2]': 0}, {'S[33]': 0, 'O[16]': -2, 'S[32]': 0, 'C[13]': 0, 'O[17]': 2, 'H[1]': 0, 'C[12]': 0, 'N[14]': 0, 'S[34]': 0, 'O[18]': 0, 'N[15]': 0, 'H[2]': 0}, {'S[33]': 0, 'O[16]': 0, 'S[32]': 0, 'C[13]': 0, 'O[17]': 0, 'H[1]': 0, 'C[12]': 0, 'N[14]': -1, 'S[34]': 0, 'O[18]': 0, 'N[15]': 1, 'H[2]': 0}, {'S[33]': 1, 'O[16]': 0, 'S[32]': -1, 'C[13]': 0, 'O[17]': 0, 'H[1]': 0, 'C[12]': 0, 'N[14]': -1, 'S[34]': 0, 'O[18]': 0, 'N[15]': 1, 'H[2]': 0}, {'S[33]': 0, 'O[16]': -1, 'S[32]': 0, 'C[13]': 0, 'O[17]': 1, 'H[1]': 0, 'C[12]': 0, 'N[14]': -1, 'S[34]': 0, 'O[18]': 0, 'N[15]': 1, 'H[2]': 0}, {'S[33]': 0, 'O[16]': 0, 'S[32]': 0, 'C[13]': 0, 'O[17]': 0, 'H[1]': 0, 'C[12]': 0, 'N[14]': -2, 'S[34]': 0, 'O[18]': 0, 'N[15]': 2, 'H[2]': 0}, {'S[33]': 0, 'O[16]': 0, 'S[32]': 0, 'C[13]': 0, 'O[17]': 0, 'H[1]': -1, 'C[12]': 0, 'N[14]': 0, 'S[34]': 0, 'O[18]': 0, 'N[15]': 0, 'H[2]': 1}, {'S[33]': 1, 'O[16]': 0, 'S[32]': -1, 'C[13]': 0, 'O[17]': 0, 'H[1]': -1, 'C[12]': 0, 'N[14]': 0, 'S[34]': 0, 'O[18]': 0, 'N[15]': 0, 'H[2]': 1}, {'S[33]': 0, 'O[16]': -1, 'S[32]': 0, 'C[13]': 0, 'O[17]': 1, 'H[1]': -1, 'C[12]': 0, 'N[14]': 0, 'S[34]': 0, 'O[18]': 0, 'N[15]': 0, 'H[2]': 1}, {'S[33]': 0, 'O[16]': 0, 'S[32]': 0, 'C[13]': 0, 'O[17]': 0, 'H[1]': -1, 'C[12]': 0, 'N[14]': -1, 'S[34]': 0, 'O[18]': 0, 'N[15]': 1, 'H[2]': 1}, {'S[33]': 0, 'O[16]': 0, 'S[32]': 0, 'C[13]': 0, 'O[17]': 0, 'H[1]': -2, 'C[12]': 0, 'N[14]': 0, 'S[34]': 0, 'O[18]': 0, 'N[15]': 0, 'H[2]': 2}, {'S[33]': 0, 'O[16]': 0, 'S[32]': 0, 'C[13]': 1, 'O[17]': 0, 'H[1]': 0, 'C[12]': -1, 'N[14]': 0, 'S[34]': 0, 'O[18]': 0, 'N[15]': 0, 'H[2]': 0}, {'S[33]': 1, 'O[16]': 0, 'S[32]': -1, 'C[13]': 1, 'O[17]': 0, 'H[1]': 0, 'C[12]': -1, 'N[14]': 0, 'S[34]': 0, 'O[18]': 0, 'N[15]': 0, 'H[2]': 0}, {'S[33]': 0, 'O[16]': -1, 'S[32]': 0, 'C[13]': 1, 'O[17]': 1, 'H[1]': 0, 'C[12]': -1, 'N[14]': 0, 'S[34]': 0, 'O[18]': 0, 'N[15]': 0, 'H[2]': 0}, {'S[33]': 0, 'O[16]': 0, 'S[32]': 0, 'C[13]': 1, 'O[17]': 0, 'H[1]': 0, 'C[12]': -1, 'N[14]': -1, 'S[34]': 0, 'O[18]': 0, 'N[15]': 1, 'H[2]': 0}, {'S[33]': 0, 'O[16]': 0, 'S[32]': 0, 'C[13]': 1, 'O[17]': 0, 'H[1]': -1, 'C[12]': -1, 'N[14]': 0, 'S[34]': 0, 'O[18]': 0, 'N[15]': 0, 'H[2]': 1}, {'S[33]': 0, 'O[16]': 0, 'S[32]': 0, 'C[13]': 2, 'O[17]': 0, 'H[1]': 0, 'C[12]': -2, 'N[14]': 0, 'S[34]': 0, 'O[18]': 0, 'N[15]': 0, 'H[2]': 0}]
    >>> all([str(LI[i][key])[:10]==str(ref[i][key])[:10] for i in range(len(LI)) for key in LI[i]])
    True
    >>> compLib['FVYHLSDLCK+++'].overlappingPeaks=={'y10 l2+': 'y10 h0+', 'y4 l1++': None, 'y2 h0+': 'y2 l2+', 'y2 h2+': None, 'y6 h2+': None, 'y3 l0++': None, 'y4 l0++': None, 'y5 h2+': None, 'y7 l1++': None, 'y9 l0+': None, 'y2 l1+': None, 'y8 h1+': None, 'y7 l0++': None, 'y10 l0+': None, 'y2 l2++': 'y2 h0++', 'y7 l0+': None, 'y5 h1+': None, 'y2 l2+': 'y2 h0+', 'y8 l1+': None, 'y5 l1+': None, 'y4 h0++': 'y4 l2++', 'y5 l1++': None, 'y4 l1+': None, 'y8 h0++': 'y8 l2++', 'y5 h0+': 'y5 l2+', 'y2 h1++': None, 'y2 l0++': None, 'y6 h0++': 'y6 l2++', 'y8 h0+': 'y8 l2+', 'y3 l1++': None, 'y8 l0+': None, 'y9 h1++': None, 'l2': 'h0', 'y10 h1++': None, 'y6 l1++': None, 'y10 h0++': 'y10 l2++', 'y8 l1++': None, 'y5 h0++': 'y5 l2++', 'y7 l2+': 'y7 h0+', 'y5 l0++': None, 'y2 l0+': None, 'y2 l1++': None, 'y7 l2++': 'y7 h0++', 'y4 h1++': None, 'y6 l1+': None, 'y4 h2++': None, 'y6 h1++': None, 'y10 l1++': None, 'y3 h1+': None, 'y4 h1+': None, 'y5 l2+': 'y5 h0+', 'y6 l2++': 'y6 h0++', 'y6 h1+': None, 'y9 h0++': 'y9 l2++', 'y9 h2+': None, 'y8 h2++': None, 'y4 h2+': None, 'h1': None, 'y3 h1++': None, 'y4 l0+': None, 'y6 l2+': 'y6 h0+', 'y4 h0+': 'y4 l2+', 'y9 h2++': None, 'y9 l1+': None, 'y9 l2++': 'y9 h0++', 'y5 h1++': None, 'y6 l0++': None, 'y7 h2++': None, 'y9 h1+': None, 'y8 l0++': None, 'y7 h0+': 'y7 l2+', 'y8 l2++': 'y8 h0++', 'y6 h2++': None, 'y10 l2++': 'y10 h0++', 'h0': 'l2', 'y3 h2++': None, 'y2 h0++': 'y2 l2++', 'y7 h0++': 'y7 l2++', 'y10 h2++': None, 'y5 h2++': None, 'y9 l1++': None, 'y6 l0+': None, 'y7 h1++': None, 'y3 h0++': 'y3 l2++', 'y8 l2+': 'y8 h0+', 'y3 l0+': None, 'y9 l0++': None, 'y7 h1+': None, 'y5 l2++': 'y5 h0++', 'y7 l1+': None, 'y3 h0+': 'y3 l2+', 'y7 h2+': None, 'y2 h2++': None, 'y3 l2+': 'y3 h0+', 'y8 h2+': None, 'y4 l2+': 'y4 h0+', 'y2 h1+': None, 'y10 h1+': None, 'y9 h0+': 'y9 l2+', 'y3 l1+': None, 'y10 h0+': 'y10 l2+', 'y9 l2+': 'y9 h0+', 'l1': None, 'y10 l0++': None, 'y8 h1++': None, 'y10 l1+': None, 'l0': None, 'y5 l0+': None, 'y3 h2+': None, 'y3 l2++': 'y3 h0++', 'y6 h0+': 'y6 l2+', 'y10 h2+': None, 'y4 l2++': 'y4 h0++', 'h2': None}
    True
    >>> str(compLib['FVYHLSDLCK+++'].precursorDistribution[0])[:10]
    '0.48945900'
    >>> str(compLib['FVYHLSDLCK+++'].precursorDistribution[1])[:10]
    '0.34951081'
    >>> str(compLib['FVYHLSDLCK+++'].precursorDistribution[2])[:10]
    '0.16103018'

    
    
    """
    pass

if __name__ == '__main__':
    
    import doctest
    doctest.testmod()
