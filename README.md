(1)main( ):在main中讀檔後存入動態記憶體，並且呼叫計算的函式時是以 陣列位址傳送。  
依照測資數目呼叫svdcmp()。   算完svd以後，w[]要以依照奇異值排好，u[]、v[]也要行列互換。   (2).myround( ):另寫一個函數處理四捨五入  
(3)svdcmp():目標是把原本的矩陣A分解成A=UWVT(旋轉VT 、伸縮W、 再旋轉U)。 U、V是正交矩陣，W是對角矩陣(sigma)，初始化wuv都是0。 判斷mn誰大誰小，如果m<n，交換mn，如果m>n，對a轉置。 把a的m行m列設給u[]，再用Jacobi奇異值分解。   u前n行裡，每行的平方和設給w，再把v算成單位矩陣，循環計算。   把w設為u中前n行，每行元素平方和的根號，然後再由大到小排好。   u和v跟w一樣的排法，u的前n行直接算 ，後面的m-n行則用Q R分解求出。   p.s.output.txt我沒有寫輸出的時候要清空檔案，資料會接在頁尾輸出
