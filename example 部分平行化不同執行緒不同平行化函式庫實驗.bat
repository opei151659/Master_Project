@echo off
set str=pla/_list.txt
for /f "tokens=*" %%a in (%str%) do (
	.\x64\Debug\Espresso-BM.exe pla\%%a -IP -IPA -ROW:1 -COL:0 -T:1 -SIG:7 -MCQE:15 -OMP >> ..\PPBM_R1_C0_T1_OMP.txt
	.\x64\Debug\Espresso-BM.exe pla\%%a -IP -IPA -ROW:1 -COL:0 -T:1 -SIG:7 -MCQE:15 -CPP >> ..\PPBM_R1_C0_T1_CPP.txt
	.\x64\Debug\Espresso-BM.exe pla\%%a -IP -IPA -ROW:1 -COL:0 -T:1 -SIG:7 -MCQE:15 -TBB >> ..\PPBM_R1_C0_T1_TBB.txt
	
	.\x64\Debug\Espresso-BM.exe pla\%%a -IP -IPA -ROW:1 -COL:0 -T:2 -SIG:7 -MCQE:15 -OMP >> ..\PPBM_R1_C0_T2_OMP.txt
	.\x64\Debug\Espresso-BM.exe pla\%%a -IP -IPA -ROW:1 -COL:0 -T:2 -SIG:7 -MCQE:15 -CPP >> ..\PPBM_R1_C0_T2_CPP.txt
	.\x64\Debug\Espresso-BM.exe pla\%%a -IP -IPA -ROW:1 -COL:0 -T:2 -SIG:7 -MCQE:15 -TBB >> ..\PPBM_R1_C0_T2_TBB.txt
	
	.\x64\Debug\Espresso-BM.exe pla\%%a -IP -IPA -ROW:1 -COL:0 -T:4 -SIG:7 -MCQE:15 -OMP >> ..\PPBM_R1_C0_T4_OMP.txt
	.\x64\Debug\Espresso-BM.exe pla\%%a -IP -IPA -ROW:1 -COL:0 -T:4 -SIG:7 -MCQE:15 -CPP >> ..\PPBM_R1_C0_T4_CPP.txt
	.\x64\Debug\Espresso-BM.exe pla\%%a -IP -IPA -ROW:1 -COL:0 -T:4 -SIG:7 -MCQE:15 -TBB >> ..\PPBM_R1_C0_T4_TBB.txt
	
	.\x64\Debug\Espresso-BM.exe pla\%%a -IP -IPA -ROW:1 -COL:0 -T:8 -SIG:7 -MCQE:15 -OMP >> ..\PPBM_R1_C0_T8_OMP.txt
	.\x64\Debug\Espresso-BM.exe pla\%%a -IP -IPA -ROW:1 -COL:0 -T:8 -SIG:7 -MCQE:15 -CPP >> ..\PPBM_R1_C0_T8_CPP.txt
	.\x64\Debug\Espresso-BM.exe pla\%%a -IP -IPA -ROW:1 -COL:0 -T:8 -SIG:7 -MCQE:15 -TBB >> ..\PPBM_R1_C0_T8_TBB.txt
	
	.\x64\Debug\Espresso-BM.exe pla\%%a -IP -IPA -ROW:1 -COL:0 -T:16 -SIG:7 -MCQE:15 -OMP >> ..\PPBM_R1_C0_T16_OMP.txt
	.\x64\Debug\Espresso-BM.exe pla\%%a -IP -IPA -ROW:1 -COL:0 -T:16 -SIG:7 -MCQE:15 -CPP >> ..\PPBM_R1_C0_T16_CPP.txt
	.\x64\Debug\Espresso-BM.exe pla\%%a -IP -IPA -ROW:1 -COL:0 -T:16 -SIG:7 -MCQE:15 -TBB >> ..\PPBM_R1_C0_T16_TBB.txt
	
	.\x64\Debug\Espresso-BM.exe pla\%%a -IP -IPA -ROW:1 -COL:0 -T:24 -SIG:7 -MCQE:15 -OMP >> ..\PPBM_R1_C0_T24_OMP.txt
	.\x64\Debug\Espresso-BM.exe pla\%%a -IP -IPA -ROW:1 -COL:0 -T:24 -SIG:7 -MCQE:15 -CPP >> ..\PPBM_R1_C0_T24_CPP.txt
	.\x64\Debug\Espresso-BM.exe pla\%%a -IP -IPA -ROW:1 -COL:0 -T:24 -SIG:7 -MCQE:15 -TBB >> ..\PPBM_R1_C0_T24_TBB.txt
)
pause
