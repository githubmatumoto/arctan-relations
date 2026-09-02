\\ -*- sh -*-
/*
  作者: 松元隆二 (matsumoto(AT)tech-i.kyutech.ac.jp)
  last updated : 2016/11/03 (JST)

  ファイル「atan_solver_check.gp」のサンプルで, 特別な式を掲載したものです。

  2015/1/13: 初版公開(maybe)
  2016/11/3: 既知の4公式は手作業でsort時に順番間違ったようなので修正。
*/

print("# 含まれるa(1/x)のxが大きいbase5 (Joerg Arndtさんの評価値)");

print("# 既知の3項式");
atan_solver([18,57,239]);
atan_solver([10,239,515]);
atan_solver([8,57,239]);
atan_solver([8,18,57]);
atan_solver([8,18,239]);

print("# 既知の4項式");
/*
作業メモ:
gzip -dc atan_table_k4.txt.gz | awk -F'[,\\[\\] ]' '{if($2>30) print $2,$3,$4,$5}' |sort -n -r | head -10
  -> 手作業でsort
*/
atan_solver([57,239,682,12943]);
atan_solver([57,68,117,239]);
atan_solver([53,57,239,4443]);
atan_solver([49,57,239,110443]);
atan_solver([43,57,117,239]);
\\atan_solver([43,57,68,239]);

print("# 既知の5項式");
/*
作業メモ:
gzip -dc atan_table_k5.txt.gz | awk -F'[,\\[\\] ]' '{if($2>100) print $2,$3,$4,$5,$6}' |sort -n -r | head -10
  -> 手作業でsort
*/
atan_solver([192,239,515,1068,173932]);
atan_solver([172,239,682,5357,12943]);
atan_solver([114,239,682,12943,740943]);
atan_solver([114,239,268,247057,740943]);
atan_solver([111,239,515,682,12943]);

print("# 既知の6項式");
/*
作業メモ:
gzip -dc atan_table_k6.txt.gz | awk -F'[,\\[\\] ]' '{if($2>260) print $2,$3,$4,$5,$6,$7}' | sort -n -r | head -10
  -> 手作業でsort
*/
atan_solver([577,682,1393,12943,32807,1049433]);
atan_solver([408,682,1393,12943,32807,1049433]);
atan_solver([408,577,682,12943,32807,1049433]);
atan_solver([355,515,682,15140,45807,89193]);
atan_solver([355,515,682,12943,45807,89193]);

print("\n# K*(pi/4) = atan(1/X) .. でKが大きい式");

/*
作業メモ:
注意:pluto-mearge-6*.DAT, conv_atan_check.plは非公開ファイル
cat pluto-mearge-6*.DAT | awk -F, '{if($3 >= 36) print $0}' | sort -n -t, -k3 -r > tmp.txt
cat tmp.txt | perl ~/2n1/search/conv-atan_check.pl -s | awk '{print "atan_solver("$1");"}' > tmp2.txt
*/
atan_solver([2,239,4193,4246,39307,390112]);
atan_solver([2,239,2855,58898,110443,4841182]);
atan_solver([2,41,463,4193,39307,390112]);
atan_solver([2,463,4193,4246,39307,390112]);
atan_solver([2,70,4193,4246,39307,390112]);
atan_solver([3,239,4193,4246,39307,390112]);
atan_solver([2,75,4193,4246,39307,390112]);

print("\n# 類似した項式 (偶然目に止まったもの)");

atan_solver([106,4443,11343,110443,595667,4841182]);
atan_solver([107,4443,11343,110443,595667,4841182]);

/* EOF */
