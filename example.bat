@echo off
set str=pla/_list.txt
for /f "tokens=*" %%a in (%str%) do (
	.\x64\Debug\Espresso-BM.exe pla\%%a pla\%%a 0 1 1 0 0 0 1 0 1 100 2 3 15 1 >> ..\Z_data_save\0910_ALL_PBM_r1_c0_T1.txt
	.\x64\Debug\Espresso-BM.exe pla\%%a pla\%%a 0 1 1 0 0 0 1 0 2 100 2 3 15 1 >> ..\Z_data_save\0910_ALL_PBM_r1_c0_T2.txt
	.\x64\Debug\Espresso-BM.exe pla\%%a pla\%%a 0 1 1 0 0 0 1 0 4 100 2 3 15 1 >> ..\Z_data_save\0910_ALL_PBM_r1_c0_T4.txt
	.\x64\Debug\Espresso-BM.exe pla\%%a pla\%%a 0 1 1 0 0 0 1 0 16 100 2 3 15 1 >> ..\Z_data_save\0910_ALL_PBM_r1_c0_T16.txt
	.\x64\Debug\Espresso-BM.exe pla\%%a pla\%%a 0 1 1 0 0 0 1 0 24 100 2 3 15 1 >> ..\Z_data_save\0910_ALL_PBM_r1_c0_T24.txt

	.\x64\Debug\Espresso-BM.exe pla\%%a pla\%%a 0 1 1 0 0 0 0 1 1 100 2 3 15 1 >> ..\Z_data_save\0910_ALL_PBM_r0_c1_T1.txt
	.\x64\Debug\Espresso-BM.exe pla\%%a pla\%%a 0 1 1 0 0 0 0 1 2 100 2 3 15 1 >> ..\Z_data_save\0910_ALL_PBM_r0_c1_T2.txt
	.\x64\Debug\Espresso-BM.exe pla\%%a pla\%%a 0 1 1 0 0 0 0 1 4 100 2 3 15 1 >> ..\Z_data_save\0910_ALL_PBM_r0_c1_T4.txt
	.\x64\Debug\Espresso-BM.exe pla\%%a pla\%%a 0 1 1 0 0 0 0 1 16 100 2 3 15 1 >> ..\Z_data_save\0910_ALL_PBM_r0_c1_T16.txt
	.\x64\Debug\Espresso-BM.exe pla\%%a pla\%%a 0 1 1 0 0 0 0 1 24 100 2 3 15 1 >> ..\Z_data_save\0910_ALL_PBM_r0_c1_T24.txt

	.\x64\Debug\Espresso-BM.exe pla\%%a pla\%%a 0 1 1 0 0 0 4 4 1 100 2 3 15 1 >> ..\Z_data_save\0910_ALL_PBM_r4_c4_T1.txt
	.\x64\Debug\Espresso-BM.exe pla\%%a pla\%%a 0 1 1 0 0 0 4 4 2 100 2 3 15 1 >> ..\Z_data_save\0910_ALL_PBM_r4_c4_T2.txt
	.\x64\Debug\Espresso-BM.exe pla\%%a pla\%%a 0 1 1 0 0 0 4 4 4 100 2 3 15 1 >> ..\Z_data_save\0910_ALL_PBM_r4_c4_T4.txt
	.\x64\Debug\Espresso-BM.exe pla\%%a pla\%%a 0 1 1 0 0 0 4 4 16 100 2 3 15 1 >> ..\Z_data_save\0910_ALL_PBM_r4_c4_T16.txt
	.\x64\Debug\Espresso-BM.exe pla\%%a pla\%%a 0 1 1 0 0 0 4 4 24 100 2 3 15 1 >> ..\Z_data_save\0910_ALL_PBM_r4_c4_T24.txt

	.\x64\Debug\Espresso-BM.exe pla\%%a pla\%%a 0 1 1 0 0 0 6 4 1 100 2 3 15 1 >> ..\Z_data_save\0910_ALL_PBM_r6_c4_T1.txt
	.\x64\Debug\Espresso-BM.exe pla\%%a pla\%%a 0 1 1 0 0 0 6 4 2 100 2 3 15 1 >> ..\Z_data_save\0910_ALL_PBM_r6_c4_T2.txt
	.\x64\Debug\Espresso-BM.exe pla\%%a pla\%%a 0 1 1 0 0 0 6 4 4 100 2 3 15 1 >> ..\Z_data_save\0910_ALL_PBM_r6_c4_T4.txt
	.\x64\Debug\Espresso-BM.exe pla\%%a pla\%%a 0 1 1 0 0 0 6 4 16 100 2 3 15 1 >> ..\Z_data_save\0910_ALL_PBM_r6_c4_T16.txt
	.\x64\Debug\Espresso-BM.exe pla\%%a pla\%%a 0 1 1 0 0 0 6 4 24 100 2 3 15 1 >> ..\Z_data_save\0910_ALL_PBM_r6_c4_T24.txt

	.\x64\Debug\Espresso-BM.exe pla\%%a pla\%%a 0 1 1 0 0 0 4 6 1 100 2 3 15 1 >> ..\Z_data_save\0910_ALL_PBM_r4_c6_T1.txt
	.\x64\Debug\Espresso-BM.exe pla\%%a pla\%%a 0 1 1 0 0 0 4 6 2 100 2 3 15 1 >> ..\Z_data_save\0910_ALL_PBM_r4_c6_T2.txt
	.\x64\Debug\Espresso-BM.exe pla\%%a pla\%%a 0 1 1 0 0 0 4 6 4 100 2 3 15 1 >> ..\Z_data_save\0910_ALL_PBM_r4_c6_T4.txt
	.\x64\Debug\Espresso-BM.exe pla\%%a pla\%%a 0 1 1 0 0 0 4 6 16 100 2 3 15 1 >> ..\Z_data_save\0910_ALL_PBM_r4_c6_T16.txt
	.\x64\Debug\Espresso-BM.exe pla\%%a pla\%%a 0 1 1 0 0 0 4 6 24 100 2 3 15 1 >> ..\Z_data_save\0910_ALL_PBM_r4_c6_T24.txt


	.\x64\Debug\Espresso-BM.exe pla\%%a pla\%%a 0 1 1 0 0 0 10 10 1 100 2 3 15 1 >> ..\Z_data_save\0910_ALL_PBM_r10_c10_T1.txt
	.\x64\Debug\Espresso-BM.exe pla\%%a pla\%%a 0 1 1 0 0 0 10 10 2 100 2 3 15 1 >> ..\Z_data_save\0910_ALL_PBM_r10_c10_T2.txt
	.\x64\Debug\Espresso-BM.exe pla\%%a pla\%%a 0 1 1 0 0 0 10 10 4 100 2 3 15 1 >> ..\Z_data_save\0910_ALL_PBM_r10_c10_T4.txt
	.\x64\Debug\Espresso-BM.exe pla\%%a pla\%%a 0 1 1 0 0 0 10 10 16 100 2 3 15 1 >> ..\Z_data_save\0910_ALL_PBM_r10_c10_T16.txt
	.\x64\Debug\Espresso-BM.exe pla\%%a pla\%%a 0 1 1 0 0 0 10 10 24 100 2 3 15 1 >> ..\Z_data_save\0910_ALL_PBM_r10_c10_T24.txt

	rem >> ..\Z_data_save\0611_OAND_ALL.txt
	rem >> ..\Z_data_save\0315pre_process_supercube.txt
)
pause
