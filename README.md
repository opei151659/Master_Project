# Master_Project
Parallel Boolean Matching based on Espresso  
  
# 安裝步驟
1.下載整個專案 並解壓縮 
2.進入external資料夾，將boost、Espresso與openapi解壓縮，解壓縮後資料夾中會有三個資料夾，如下:  
  (如無法順利解壓縮請到下面使用library下載並放到對應位置)
  external  
    ---- boost_1_82_0  
    ---- Espresso  
    ---- openapi-tbb-2021.10.0  
3. 安裝visual studio 2022 後，點Espresso-BM-ver6.0.sin檔 開啟專案  
4. 開啟專案後即可執行  
  
  
# 測試電路
 pla 資料夾為各式pla電路，內有兩個完整的benchmark壓縮檔 MCNC與ISCAS 

  
# 輸入參數說明:
name.pla -IP -IPA

如需比對兩個電路: (必須人工確定兩個電路的輸入與輸出個數相同)
name.pla name.pla -IP -IPA

-h -H -help -HELP 會出現此文檔

NUM代表整數 FNUM代表小數  
-IN:NUM   擴展輸入的input個數，數字小於原本大小包持原樣，最大為128 (0 ~ 128)   
-IP     考慮input permutation  
-IPA    考慮input phase assignment  
-OP     考慮output permutation  
-OPA    考慮output phase assignment  
-RIP    打亂輸入random input permutation  
-RIPA   打亂輸入random input phase assignment  
-ROP    打亂輸入random output permutation  
-ROPA   打亂輸出random output phase assignment  
-RM:NUM   移除on-set 與off-set 的百分比 (0.0 ~ 1.1) ex: 0.1表示移除10%  
-ROW:NUM  row(on-set) 要分割的份數 0:表示分割成#row 份 1: 表示無分割  
-COL:NUM  column(off-set) 要分割的份數 0:表示分割成#col 份 1: 表示無分割  
-T:NUM    使用執行緒的最大個數 (最小為1) 最好<=CPU的最大值行緒個數有做限制  
-CHK:NUM    執行相同任務時每個執行緒所分配的區塊大小 (最小為1) 100適合絕大多數情況  
-OMP    使用的平行化函式庫OMP 不可與'-CPP'與'-TBB'同時使用  
-CPP      使用的平行化函式庫CPP 不可與'-OMP'與'-TBB'同時使用  
-TBB      使用的平行化函式庫TBB 不可與'-OMP'與'-CPP'同時使用  
-SIG:NUM  選擇使用的輸入特徵值(輸入為10進位) 表示為 0b0001: check_unate, 0b0010 = check_ENE, 0b0100 = check_KSIG, 0b1000 = check_cofactor 可以同時使用多個 (15(  1111))  
-MCQE:NUM 選擇使用的MCQE RULE(輸入為10進位) 表示為 0b0001: Rule1, 0b0010 = Rule2, 0b0100 = Rule3, 0b1000 = Rule4 可以同時使用多個 (15(1111))  
-NTP    不使用執行緒池進行平行化，只可與'-CPP'以起始用    
-ALLP   使用全部平行化，未使用則為部分平行化，只可與'-CPP'以起始用  

### 以下4個指令會自動將設定對應實驗參數 僅需使用以下單個指令 建議不要與其他指令共用       
-OutputSIG  執行輸出特徵值實驗  (-IP -IPA -OP -T:1)
-InputSIG 執行輸入特徵值實驗  (-IP -IPA -T:1)
-ALLMCQE  執行MCQE所有組合實驗  (-IP -IPA -T:1)
-WINDOWSIZE 執行複數積項配對window size = 2實驗  (-IP -IPA -T:1)
  
預設參數  
-CPP -SIG:7 -MCQE:15 -T:1 -ROW:1 -COL:1 -RM:0  
  
# 範例  
.\x64\Debug\Espresso-BM.exe pla\alu2.pla  -IP -IPA -T:24   

# 資料庫連結設定 (下載完整專案則不須手動設定)
1.在vs屬性頁 -> C++ -> 一般 -> 其他Include目錄 中加入  
Espresso: Espresso\src; 的路徑  
OneAPI: oneapi-tbb-2021.10.0\include; 的路徑  
Boost: boost_1_82_0; 的路徑
  
2.在vs屬性頁 -> 連結器 -> 一般 -> 其他程式庫目錄 中加入  
Espresso: Espresso\Debug; 的路徑  
OneAPI: oneapi-tbb-2021.10.0\lib\intel64\vc14; 的路徑  
Boost: boost_1_82_0\lib; 的路徑
  
3.在vs屬性頁 -> 連結器 -> 輸入 -> 其他相依性中加入 ESPRESSO.lib;  
  
4.在vs屬性頁 -> C++ -> 語言 勾選開啟MP支援 啟用OpenMP  
5.在vs屬性頁 -> C++ -> 語言 C++ 語言表準改成C++ 17 or C++ 20                       
  
# 可能的Debug 排除  (下載完整專案則不須手動設定) 
1. 在vs屬性頁 -> C++ -> 語言 -> 一致性模式 改成否(/permissive)  
2. 在vs屬性頁 -> C++ -> 前置處理器 -> 前置處理器定義加入 _CRT_SECURE_NO_WARNINGS;  

# 使用library
## Espresso
https://github.com/Gigantua/Espresso)https://github.com/Gigantua/Espresso
## boost
 https://www.boost.org/ 使用的版本為boost_1_82_0
## IntelTBB
https://github.com/oneapi-src/oneTBB 使用的版本為oneapi-tbb-2021.9.0
