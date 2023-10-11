# Mater_Project
Parallel Boolean Matching based on Espresso  
  
# 安裝步驟
1.下載整個專案 並解壓縮  
2.進入external資料夾，將boost、Espresso與openapi解壓縮，解壓縮後資料夾中會有三個資料夾，如下:  
  external  
    ---- boost_1_82_0  
    ---- Espresso  
    ---- openapi-tbb-2021.10.0  
3. 安裝visual studio 2022 後，點Espresso-BM-ver6.0.sin檔 開啟專案  
4. 開啟專案後即可執行  
  
  
# 測試電路
 pla 資料夾為各式pla電路，內有兩個完整的benchmark壓縮檔 MCNC與ISCAS 

  
# 輸入參數說明:
name.pla name.pla A B C D E F G H I J K L M N O   
A: 是否要擴展輸入的input個數，數字小於原本大小包持原樣，最大為128  (0 ~ 128)  
B: 是否有input permutation (是: 1/否: 0)  
C: 是否有input phase assignment (是: 1/否: 0)  
D: 是否有output permutation (是: 1/否: 0)  
E: 是否有output phase assignment (是: 1/否: 0)  
F: 移除on-set 與off-set 的百分比  (0.0 ~ 1.1) ex: 0.1表示移除10%  
G: row(on-set) 要分割的份數 0:表示分割成#row 份 1: 表示無分割  
H: column(off-set) 要分割的份數 0:表示分割成#col 份 1: 表示無分割  
I: 使用執行緒的最大個數 (最小為1) 最好<=CPU的最大值行緒個數 否則會有額外負擔  
J: 執行相同任務時每個執行緒所分配的區塊大小 (最小為1) 100適合絕大多數情況  
K: 使用的平行化函式庫 1:OMP 2:CPP 3:TBB (僅在部分平行化時可更改) 其餘狀況需都設為2(CPP)  
L: 選擇使用的輸入特徵值  0b0001: check_unate, 0b0010 = check_ENE, 0b0100 = check_KSIG, 0b1000 = check_cofactor 可以同時使用多個 (15(1111))  
M: 選擇使用的MCQE  0b0001: Rule1, 0b0010 = Rule2, 0b0100 = Rule3, 0b1000 = Rule4 可以同時使用多個 (15(1111))  
O: 是否執行全部平行化 (是: 1/否: 0)  

# 範例
.\x64\Debug\Espresso-BM.exe pla\alu2.pla pla\alu2.pla 0 1 1 0 0 0 1 1 1 100 2 7 15 0  

