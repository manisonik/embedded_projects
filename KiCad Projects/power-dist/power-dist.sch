EESchema Schematic File Version 4
LIBS:power-dist-cache
EELAYER 26 0
EELAYER END
$Descr A4 11693 8268
encoding utf-8
Sheet 1 1
Title ""
Date ""
Rev ""
Comp ""
Comment1 ""
Comment2 ""
Comment3 ""
Comment4 ""
$EndDescr
$Comp
L power-dist-rescue:MP1584-Regulator_Switching U1
U 1 1 5C1761B2
P 2700 1900
F 0 "U1" H 2400 2400 50  0000 C CNN
F 1 "MP1584" H 2900 2400 50  0000 C CNN
F 2 "Housings_SOIC:SOIC-8-1EP_3.9x4.9mm_Pitch1.27mm" H 3800 1400 50  0001 C CNN
F 3 "https://www.monolithicpower.com/pub/media/document/MP1584_r1.0.pdf" H 4200 1300 50  0001 C CNN
	1    2700 1900
	1    0    0    -1  
$EndComp
Wire Wire Line
	2700 1250 2700 1150
Wire Wire Line
	3350 1150 3350 1550
Wire Wire Line
	3250 1550 3350 1550
Wire Wire Line
	3350 1550 3500 1550
Connection ~ 3350 1550
$Comp
L Device:C_Small C2
U 1 1 5C17631C
P 3050 1150
F 0 "C2" V 2821 1150 50  0000 C CNN
F 1 "100nF" V 2912 1150 50  0000 C CNN
F 2 "Capacitors_SMD:C_0201" H 3050 1150 50  0001 C CNN
F 3 "~" H 3050 1150 50  0001 C CNN
	1    3050 1150
	0    1    1    0   
$EndComp
Wire Wire Line
	2700 1150 2950 1150
Wire Wire Line
	3150 1150 3350 1150
Wire Wire Line
	2150 2250 2050 2250
Wire Wire Line
	2050 2250 2050 2350
$Comp
L Device:R_Small_US R3
U 1 1 5C17676B
P 2050 2450
F 0 "R3" H 2118 2496 50  0000 L CNN
F 1 "100k" H 2118 2405 50  0000 L CNN
F 2 "Resistors_SMD:R_0201" H 2050 2450 50  0001 C CNN
F 3 "~" H 2050 2450 50  0001 C CNN
	1    2050 2450
	1    0    0    -1  
$EndComp
$Comp
L power:GND #PWR0101
U 1 1 5C176877
P 2050 2650
F 0 "#PWR0101" H 2050 2400 50  0001 C CNN
F 1 "GND" H 2055 2477 50  0000 C CNN
F 2 "" H 2050 2650 50  0001 C CNN
F 3 "" H 2050 2650 50  0001 C CNN
	1    2050 2650
	1    0    0    -1  
$EndComp
Wire Wire Line
	2050 2550 2050 2650
Wire Wire Line
	2700 2550 2700 2650
$Comp
L power:GND #PWR0102
U 1 1 5C176ADA
P 2700 2650
F 0 "#PWR0102" H 2700 2400 50  0001 C CNN
F 1 "GND" H 2705 2477 50  0000 C CNN
F 2 "" H 2700 2650 50  0001 C CNN
F 3 "" H 2700 2650 50  0001 C CNN
	1    2700 2650
	1    0    0    -1  
$EndComp
Wire Wire Line
	2700 2550 2900 2550
Wire Wire Line
	2900 2550 2900 2650
Connection ~ 2700 2550
$Comp
L power:GNDREF #PWR0103
U 1 1 5C176E85
P 2900 2650
F 0 "#PWR0103" H 2900 2400 50  0001 C CNN
F 1 "GNDREF" H 3100 2650 50  0000 C CNN
F 2 "" H 2900 2650 50  0001 C CNN
F 3 "" H 2900 2650 50  0001 C CNN
	1    2900 2650
	1    0    0    -1  
$EndComp
$Comp
L Device:C_Small C3
U 1 1 5C177FF9
P 3550 2500
F 0 "C3" H 3642 2546 50  0000 L CNN
F 1 "150pF" H 3642 2455 50  0000 L CNN
F 2 "Capacitors_SMD:C_0201" H 3550 2500 50  0001 C CNN
F 3 "~" H 3550 2500 50  0001 C CNN
	1    3550 2500
	1    0    0    -1  
$EndComp
Wire Wire Line
	4000 2250 4000 2450
$Comp
L Device:C_Small C4
U 1 1 5C1783E7
P 4000 2550
F 0 "C4" H 4092 2596 50  0000 L CNN
F 1 "NS" H 4092 2505 50  0000 L CNN
F 2 "Capacitors_SMD:C_0201" H 4000 2550 50  0001 C CNN
F 3 "~" H 4000 2550 50  0001 C CNN
	1    4000 2550
	1    0    0    -1  
$EndComp
$Comp
L Device:R_Small_US R5
U 1 1 5C178643
P 3550 2750
F 0 "R5" H 3618 2796 50  0000 L CNN
F 1 "100k" H 3618 2705 50  0000 L CNN
F 2 "Resistors_SMD:R_0201" H 3550 2750 50  0001 C CNN
F 3 "~" H 3550 2750 50  0001 C CNN
	1    3550 2750
	1    0    0    -1  
$EndComp
$Comp
L power:GND #PWR0106
U 1 1 5C1786DD
P 3550 2850
F 0 "#PWR0106" H 3550 2600 50  0001 C CNN
F 1 "GND" H 3555 2677 50  0000 C CNN
F 2 "" H 3550 2850 50  0001 C CNN
F 3 "" H 3550 2850 50  0001 C CNN
	1    3550 2850
	1    0    0    -1  
$EndComp
Wire Wire Line
	3550 2850 4000 2850
Wire Wire Line
	4000 2850 4000 2650
Connection ~ 3550 2850
$Comp
L pspice:INDUCTOR L1
U 1 1 5C178B29
P 4000 1550
F 0 "L1" H 4000 1765 50  0000 C CNN
F 1 "15uH" H 4000 1674 50  0000 C CNN
F 2 "Inductors_SMD:L_1210" H 4000 1550 50  0001 C CNN
F 3 "~" H 4000 1550 50  0001 C CNN
	1    4000 1550
	1    0    0    -1  
$EndComp
$Comp
L Device:D_Zener_Small D1
U 1 1 5C178C3A
P 3500 1650
F 0 "D1" V 3454 1718 50  0000 L CNN
F 1 "Zener" V 3545 1718 50  0000 L CNN
F 2 "Diodes_SMD:D_0603" V 3500 1650 50  0001 C CNN
F 3 "~" V 3500 1650 50  0001 C CNN
	1    3500 1650
	0    1    1    0   
$EndComp
Connection ~ 3500 1550
Wire Wire Line
	3500 1550 3750 1550
Wire Wire Line
	3250 1900 3500 1900
$Comp
L Device:R_Small_US R4
U 1 1 5C179327
P 3500 2000
F 0 "R4" H 3568 2046 50  0000 L CNN
F 1 "40.2k" H 3568 1955 50  0000 L CNN
F 2 "Resistors_SMD:R_0201" H 3500 2000 50  0001 C CNN
F 3 "~" H 3500 2000 50  0001 C CNN
	1    3500 2000
	1    0    0    -1  
$EndComp
Wire Wire Line
	4250 1550 4400 1550
$Comp
L Device:C_Small C5
U 1 1 5C179B11
P 4400 1650
F 0 "C5" H 4492 1696 50  0000 L CNN
F 1 "22uF" H 4492 1605 50  0000 L CNN
F 2 "Capacitors_SMD:C_0201" H 4400 1650 50  0001 C CNN
F 3 "~" H 4400 1650 50  0001 C CNN
	1    4400 1650
	1    0    0    -1  
$EndComp
$Comp
L power:GNDREF #PWR0107
U 1 1 5C179BA1
P 4400 1750
F 0 "#PWR0107" H 4400 1500 50  0001 C CNN
F 1 "GNDREF" H 4200 1750 50  0000 C CNN
F 2 "" H 4400 1750 50  0001 C CNN
F 3 "" H 4400 1750 50  0001 C CNN
	1    4400 1750
	1    0    0    -1  
$EndComp
$Comp
L power:GNDREF #PWR0108
U 1 1 5C179C35
P 3500 1750
F 0 "#PWR0108" H 3500 1500 50  0001 C CNN
F 1 "GNDREF" H 3300 1750 50  0000 C CNN
F 2 "" H 3500 1750 50  0001 C CNN
F 3 "" H 3500 1750 50  0001 C CNN
	1    3500 1750
	1    0    0    -1  
$EndComp
$Comp
L power:GNDREF #PWR0109
U 1 1 5C179C6F
P 3500 2100
F 0 "#PWR0109" H 3500 1850 50  0001 C CNN
F 1 "GNDREF" H 3300 2100 50  0000 C CNN
F 2 "" H 3500 2100 50  0001 C CNN
F 3 "" H 3500 2100 50  0001 C CNN
	1    3500 2100
	1    0    0    -1  
$EndComp
Wire Wire Line
	4400 1550 4800 1550
Connection ~ 4400 1550
Wire Wire Line
	4800 1900 4800 1550
Connection ~ 3500 1900
$Comp
L Device:R_Small_US R6
U 1 1 5C17AD50
P 4150 1900
F 0 "R6" V 4250 1900 50  0000 C CNN
F 1 "210k" V 4350 1900 50  0000 C CNN
F 2 "Resistors_SMD:R_0201" H 4150 1900 50  0001 C CNN
F 3 "~" H 4150 1900 50  0001 C CNN
	1    4150 1900
	0    1    1    0   
$EndComp
Wire Wire Line
	4250 1900 4800 1900
Wire Wire Line
	3500 1900 4050 1900
Wire Wire Line
	4800 1550 5000 1550
Connection ~ 4800 1550
Text GLabel 5000 1550 2    50   Output ~ 0
5V_OUT
Wire Wire Line
	3250 2250 3550 2250
Wire Wire Line
	3550 2250 3550 2400
Connection ~ 3550 2250
Wire Wire Line
	3550 2250 4000 2250
Wire Wire Line
	3550 2600 3550 2650
$Comp
L power:GNDREF #PWR0104
U 1 1 5C18BFCD
P 2150 2050
F 0 "#PWR0104" H 2150 1800 50  0001 C CNN
F 1 "GNDREF" V 2250 2200 50  0000 R CNN
F 2 "" H 2150 2050 50  0001 C CNN
F 3 "" H 2150 2050 50  0001 C CNN
	1    2150 2050
	0    1    1    0   
$EndComp
$Comp
L power-dist-rescue:MP1584-Regulator_Switching U2
U 1 1 5C181E0C
P 2550 4200
F 0 "U2" H 2250 4700 50  0000 C CNN
F 1 "MP1584" H 2750 4700 50  0000 C CNN
F 2 "Housings_SOIC:SOIC-8-1EP_3.9x4.9mm_Pitch1.27mm" H 3650 3700 50  0001 C CNN
F 3 "https://www.monolithicpower.com/pub/media/document/MP1584_r1.0.pdf" H 4050 3600 50  0001 C CNN
	1    2550 4200
	1    0    0    -1  
$EndComp
Wire Wire Line
	2550 3550 2550 3450
Wire Wire Line
	3200 3450 3200 3850
Wire Wire Line
	3100 3850 3200 3850
Wire Wire Line
	3200 3850 3350 3850
Connection ~ 3200 3850
$Comp
L Device:C_Small C1
U 1 1 5C181E18
P 2900 3450
F 0 "C1" V 2671 3450 50  0000 C CNN
F 1 "100nF" V 2762 3450 50  0000 C CNN
F 2 "Capacitors_SMD:C_0201" H 2900 3450 50  0001 C CNN
F 3 "~" H 2900 3450 50  0001 C CNN
	1    2900 3450
	0    1    1    0   
$EndComp
Wire Wire Line
	2550 3450 2800 3450
Wire Wire Line
	3000 3450 3200 3450
Wire Wire Line
	2000 4550 1900 4550
Wire Wire Line
	1900 4550 1900 4650
$Comp
L Device:R_Small_US R1
U 1 1 5C181E24
P 1900 4750
F 0 "R1" H 1968 4796 50  0000 L CNN
F 1 "200k" H 1968 4705 50  0000 L CNN
F 2 "Resistors_SMD:R_0201" H 1900 4750 50  0001 C CNN
F 3 "~" H 1900 4750 50  0001 C CNN
	1    1900 4750
	1    0    0    -1  
$EndComp
$Comp
L power:GND #PWR0105
U 1 1 5C181E2B
P 1900 4950
F 0 "#PWR0105" H 1900 4700 50  0001 C CNN
F 1 "GND" H 1905 4777 50  0000 C CNN
F 2 "" H 1900 4950 50  0001 C CNN
F 3 "" H 1900 4950 50  0001 C CNN
	1    1900 4950
	1    0    0    -1  
$EndComp
Wire Wire Line
	1900 4850 1900 4950
Wire Wire Line
	2550 4850 2550 4950
$Comp
L power:GND #PWR0110
U 1 1 5C181E33
P 2550 4950
F 0 "#PWR0110" H 2550 4700 50  0001 C CNN
F 1 "GND" H 2555 4777 50  0000 C CNN
F 2 "" H 2550 4950 50  0001 C CNN
F 3 "" H 2550 4950 50  0001 C CNN
	1    2550 4950
	1    0    0    -1  
$EndComp
Wire Wire Line
	2550 4850 2750 4850
Wire Wire Line
	2750 4850 2750 4950
Connection ~ 2550 4850
$Comp
L power:GNDREF #PWR0111
U 1 1 5C181E3C
P 2750 4950
F 0 "#PWR0111" H 2750 4700 50  0001 C CNN
F 1 "GNDREF" H 2950 4950 50  0000 C CNN
F 2 "" H 2750 4950 50  0001 C CNN
F 3 "" H 2750 4950 50  0001 C CNN
	1    2750 4950
	1    0    0    -1  
$EndComp
$Comp
L Device:C_Small C6
U 1 1 5C181E62
P 3400 4800
F 0 "C6" H 3492 4846 50  0000 L CNN
F 1 "220pF" H 3492 4755 50  0000 L CNN
F 2 "Capacitors_SMD:C_0201" H 3400 4800 50  0001 C CNN
F 3 "~" H 3400 4800 50  0001 C CNN
	1    3400 4800
	1    0    0    -1  
$EndComp
Wire Wire Line
	3850 4550 3850 4750
$Comp
L Device:C_Small C7
U 1 1 5C181E6A
P 3850 4850
F 0 "C7" H 3942 4896 50  0000 L CNN
F 1 "NS" H 3942 4805 50  0000 L CNN
F 2 "Capacitors_SMD:C_0201" H 3850 4850 50  0001 C CNN
F 3 "~" H 3850 4850 50  0001 C CNN
	1    3850 4850
	1    0    0    -1  
$EndComp
$Comp
L Device:R_Small_US R7
U 1 1 5C181E71
P 3400 5050
F 0 "R7" H 3468 5096 50  0000 L CNN
F 1 "68.1k" H 3468 5005 50  0000 L CNN
F 2 "Resistors_SMD:R_0201" H 3400 5050 50  0001 C CNN
F 3 "~" H 3400 5050 50  0001 C CNN
	1    3400 5050
	1    0    0    -1  
$EndComp
$Comp
L power:GND #PWR0112
U 1 1 5C181E78
P 3400 5150
F 0 "#PWR0112" H 3400 4900 50  0001 C CNN
F 1 "GND" H 3405 4977 50  0000 C CNN
F 2 "" H 3400 5150 50  0001 C CNN
F 3 "" H 3400 5150 50  0001 C CNN
	1    3400 5150
	1    0    0    -1  
$EndComp
Wire Wire Line
	3400 5150 3850 5150
Wire Wire Line
	3850 5150 3850 4950
Connection ~ 3400 5150
$Comp
L pspice:INDUCTOR L2
U 1 1 5C181E81
P 3850 3850
F 0 "L2" H 3850 4065 50  0000 C CNN
F 1 "10uH" H 3850 3974 50  0000 C CNN
F 2 "Inductors_SMD:L_1210" H 3850 3850 50  0001 C CNN
F 3 "~" H 3850 3850 50  0001 C CNN
	1    3850 3850
	1    0    0    -1  
$EndComp
$Comp
L Device:D_Zener_Small D2
U 1 1 5C181E88
P 3350 3950
F 0 "D2" V 3304 4018 50  0000 L CNN
F 1 "Zener" V 3395 4018 50  0000 L CNN
F 2 "Diodes_SMD:D_0603" V 3350 3950 50  0001 C CNN
F 3 "~" V 3350 3950 50  0001 C CNN
	1    3350 3950
	0    1    1    0   
$EndComp
Connection ~ 3350 3850
Wire Wire Line
	3350 3850 3600 3850
Wire Wire Line
	3100 4200 3350 4200
$Comp
L Device:R_Small_US R2
U 1 1 5C181E92
P 3350 4300
F 0 "R2" H 3418 4346 50  0000 L CNN
F 1 "40.2k" H 3418 4255 50  0000 L CNN
F 2 "Resistors_SMD:R_0201" H 3350 4300 50  0001 C CNN
F 3 "~" H 3350 4300 50  0001 C CNN
	1    3350 4300
	1    0    0    -1  
$EndComp
Wire Wire Line
	4100 3850 4250 3850
$Comp
L Device:C_Small C8
U 1 1 5C181E9A
P 4250 3950
F 0 "C8" H 4342 3996 50  0000 L CNN
F 1 "22uF" H 4342 3905 50  0000 L CNN
F 2 "Capacitors_SMD:C_0201" H 4250 3950 50  0001 C CNN
F 3 "~" H 4250 3950 50  0001 C CNN
	1    4250 3950
	1    0    0    -1  
$EndComp
$Comp
L power:GNDREF #PWR0113
U 1 1 5C181EA1
P 4250 4050
F 0 "#PWR0113" H 4250 3800 50  0001 C CNN
F 1 "GNDREF" H 4050 4050 50  0000 C CNN
F 2 "" H 4250 4050 50  0001 C CNN
F 3 "" H 4250 4050 50  0001 C CNN
	1    4250 4050
	1    0    0    -1  
$EndComp
$Comp
L power:GNDREF #PWR0114
U 1 1 5C181EA7
P 3350 4050
F 0 "#PWR0114" H 3350 3800 50  0001 C CNN
F 1 "GNDREF" H 3150 4050 50  0000 C CNN
F 2 "" H 3350 4050 50  0001 C CNN
F 3 "" H 3350 4050 50  0001 C CNN
	1    3350 4050
	1    0    0    -1  
$EndComp
$Comp
L power:GNDREF #PWR0115
U 1 1 5C181EAD
P 3350 4400
F 0 "#PWR0115" H 3350 4150 50  0001 C CNN
F 1 "GNDREF" H 3150 4400 50  0000 C CNN
F 2 "" H 3350 4400 50  0001 C CNN
F 3 "" H 3350 4400 50  0001 C CNN
	1    3350 4400
	1    0    0    -1  
$EndComp
Wire Wire Line
	4250 3850 4650 3850
Connection ~ 4250 3850
Wire Wire Line
	4650 4200 4650 3850
Connection ~ 3350 4200
$Comp
L Device:R_Small_US R8
U 1 1 5C181EB7
P 4000 4200
F 0 "R8" V 4100 4200 50  0000 C CNN
F 1 "124k" V 4200 4200 50  0000 C CNN
F 2 "Resistors_SMD:R_0201" H 4000 4200 50  0001 C CNN
F 3 "~" H 4000 4200 50  0001 C CNN
	1    4000 4200
	0    1    1    0   
$EndComp
Wire Wire Line
	4100 4200 4650 4200
Wire Wire Line
	3350 4200 3900 4200
Wire Wire Line
	4650 3850 4850 3850
Connection ~ 4650 3850
Text GLabel 4850 3850 2    50   Output ~ 0
3.3V_OUT
Wire Wire Line
	3100 4550 3400 4550
Wire Wire Line
	3400 4550 3400 4700
Connection ~ 3400 4550
Wire Wire Line
	3400 4550 3850 4550
Wire Wire Line
	3400 4900 3400 4950
$Comp
L power:GNDREF #PWR0116
U 1 1 5C181ED9
P 2000 4350
F 0 "#PWR0116" H 2000 4100 50  0001 C CNN
F 1 "GNDREF" V 2100 4500 50  0000 R CNN
F 2 "" H 2000 4350 50  0001 C CNN
F 3 "" H 2000 4350 50  0001 C CNN
	1    2000 4350
	0    1    1    0   
$EndComp
Wire Notes Line
	500  3100 6300 3100
Text Notes 550  600  0    50   ~ 0
5V
$Comp
L Device:R_Small_US R9
U 1 1 5C1901AF
P 8050 1600
F 0 "R9" H 8118 1646 50  0000 L CNN
F 1 "24.9k" H 8118 1555 50  0000 L CNN
F 2 "Resistors_SMD:R_0201" H 8050 1600 50  0001 C CNN
F 3 "~" H 8050 1600 50  0001 C CNN
	1    8050 1600
	1    0    0    -1  
$EndComp
$Comp
L power:GNDREF #PWR0117
U 1 1 5C1901B6
P 7700 1150
F 0 "#PWR0117" H 7700 900 50  0001 C CNN
F 1 "GNDREF" H 7500 1150 50  0000 C CNN
F 2 "" H 7700 1150 50  0001 C CNN
F 3 "" H 7700 1150 50  0001 C CNN
	1    7700 1150
	1    0    0    -1  
$EndComp
$Comp
L Device:C_Small C9
U 1 1 5C1901BD
P 7700 1050
F 0 "C9" H 7792 1096 50  0000 L CNN
F 1 "10uF" H 7792 1005 50  0000 L CNN
F 2 "Capacitors_SMD:CP_Elec_4x4.5" H 7700 1050 50  0001 C CNN
F 3 "~" H 7700 1050 50  0001 C CNN
	1    7700 1050
	1    0    0    -1  
$EndComp
$Comp
L power:GND #PWR0118
U 1 1 5C1901C5
P 8050 1700
F 0 "#PWR0118" H 8050 1450 50  0001 C CNN
F 1 "GND" H 8055 1527 50  0000 C CNN
F 2 "" H 8050 1700 50  0001 C CNN
F 3 "" H 8050 1700 50  0001 C CNN
	1    8050 1700
	1    0    0    -1  
$EndComp
Wire Wire Line
	8050 1050 8050 1100
Wire Wire Line
	8050 1050 8150 1050
$Comp
L Device:R_Small_US R10
U 1 1 5C1901CD
P 8050 1200
F 0 "R10" H 8118 1246 50  0000 L CNN
F 1 "100k" H 8118 1155 50  0000 L CNN
F 2 "Resistors_SMD:R_0201" H 8050 1200 50  0001 C CNN
F 3 "~" H 8050 1200 50  0001 C CNN
	1    8050 1200
	1    0    0    -1  
$EndComp
Wire Wire Line
	8050 1300 8050 1400
Wire Wire Line
	8050 1400 8150 1400
Text GLabel 8150 1050 2    50   Output ~ 0
VIN
Text GLabel 8150 1400 2    50   Output ~ 0
EN
Text GLabel 2150 1550 0    50   Input ~ 0
RAW
Text GLabel 2150 1900 0    50   Input ~ 0
EN
Text GLabel 2000 3850 0    50   Input ~ 0
RAW
Text GLabel 2000 4200 0    50   Input ~ 0
EN
Wire Notes Line
	6300 5450 500  5450
Wire Notes Line
	6300 500  6300 5450
Text Notes 550  3250 0    50   ~ 0
3.3V
Text GLabel 7000 1200 2    50   Output ~ 0
RAW
Text GLabel 7700 950  1    50   Input ~ 0
RAW
Text GLabel 8050 950  1    50   Input ~ 0
RAW
Wire Wire Line
	8050 1400 8050 1500
Connection ~ 8050 1400
Wire Wire Line
	8050 950  8050 1050
Connection ~ 8050 1050
$Comp
L Connector:Conn_01x01_Male J1
U 1 1 5C5AB8F1
P 6800 1200
F 0 "J1" H 6906 1378 50  0000 C CNN
F 1 "Conn_01x01_Male" H 6950 1100 50  0000 C CNN
F 2 "Wire_Pads:SolderWirePad_4xSquare_0-8mmDrill" H 6800 1200 50  0001 C CNN
F 3 "~" H 6800 1200 50  0001 C CNN
	1    6800 1200
	1    0    0    -1  
$EndComp
$EndSCHEMATC