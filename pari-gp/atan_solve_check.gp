\\ -*- gp-script -*-
/*
  ライセンス: GPL
  作者: 松元隆二 (matsumoto(AT)tech-i.kyutech.ac.jp)
  last updated : 2016/11/11 (JST)  

  1. arctan関係式の係数を計算するプログラム
       - atan_solver()

  2. 既知の公式のチェックをするプログラム
       - atan_check()
        
  PARI/GPが必要です。
    http://pari.math.u-bordeaux.fr/


2014/12/13, atan_check()を追加
2014/12/13, Web掲載版:rev1543
2015/1/5, atan_solver0()を追加
2015/1/5, atan_solver2()を追加
2015/1/5, collect_prime()を整数の逆数以外も可能に拡張
2015/1/5, atan_solver0()を整数の逆数以外も可能に拡張
          自作アルゴリズムを採用しているため、誤る可能性が否定できません。
          (私は整数論はよくわからない。)
2015/1/5, collect_prime()を巨大な整数に対して高速化。簡単に素因数
          分解できなかった場合は擬似素数とみなすように変更。
          この変更により、計算を誤る可能性があります。
           特に2重ベクトルの式の判断で誤判定の可能性が否定できない。
2015/1/6, atan_solver0()で、素因数分解で擬似素数を採用した場合は
          素数リストに[pseudo]を表示するように変更。
2015/1/8, 実行例を別ファイル atan_solver_check_sample.gpに分離。

2015/1/12, Web掲載版:rev1560
2015/1/14, Web掲載版:rev1568
2015/1/18, プログラムの整理
2015/1/18, 旧形式への変換で
              atan_table_k7ov.txt.gz/atan_table_u-all.txt.gz + atan_check
           の組み合わせでメモリ不足エラーになる修正。原因は複素数乗
	   算の絶対値が大きくなりすぎていたためなので、閾値設置して戦略
	   変更。
2015/1/18, Web掲載版:rev1574
2016/5/5, 一部バグ修正。vecsortを追加。
           x setsearch(XX, 0)
           o setsearch(vecsort(XX), 0)
2016/7/20, Web掲載版:rev1644
           ヘッダにemacsのgp-script-mode用の記述追加しています。
2016/11/11, 別表のarctan関係式の全リストのリンク先を修正。計算時間を修正。

(編集作業メモ: atan_table更新時はリンク先の更新と計算時間の更新を忘れずに)
*/

/*
  1. arctan関係式の係数を計算するプログラムです。

  関数atan_solver()で係数を求めます。

  例えば

    K * pi/4  = k1 * atan(1/x1) + k2 * atan(1/x2) + k3 * atan(1/x3).

  という式で、x1,x2,x3を与えると K, k1, k2, k3を算出します。

   第一引数でxのリストを渡します。

     atan_solver([xのリスト])

  既知の公式を例に説明します。ガウスの公式

     pi/4 = +12arctan(1/18)+8arctan(1/57)-5arctan(1/239).

  で x1 =18, x2 = 57, x3 = 239なので、

     atan_solver([18,57,239]);
     -> 出力: 8933, 1, +12a(1/18)+8a(1/57)-5a(1/239)., [5,13].
        (出力フォーマットは atan_table.html参照ください。)

  とすれば、算出されます。

  もし1/xが整数の逆数以外の場合は、その逆数を入れてください。

  例えば

       pi/4 = +5arctan(29/278)+7arctan(3/79). 

  の場合、29/5の逆数は278/29, 3/79の逆数は79/3なので

      atan_solver([278/29,79/3]);
      -> 出力: 8613, 1, +5a(29/278)+7a(3/79)., [5].

  立式に成功しなかった場合のエラーがいくつかあります。サンプルファイル
  「asc_sample.gp」を参照ください。

  2. 既知の公式のチェックをするプログラム

　関数atan_check()で既知の公式のチェックをします。ただし、2重ベクトルの
  式も正常と認識します。

  例えば
    K * pi/4  = k1 * atan(1/x1) + k2 * atan(1/x2) + k3 * atan(1/x3).
    
　という式があった場合、有効な式かをチェックします。

   第一引数でxのリストを渡します。 第二引数でkのリストを渡します。

   atan_solver([xのリスト],[kのリスト])

  既知の公式を例に説明します。ガウスの公式

     pi/4 = +12arctan(1/18)+8arctan(1/57)-5arctan(1/239).

  で x1 =18, x2 = 57, x3 = 239, k1=12, k2=8, k3=-5 なので、

     atan_check([18,57,239],[12,8,-5]);
     -> 出力: 8933, 1, +12a(1/18)+8a(1/57)-5a(1/239).
       (出力フォーマットは atan_table.html参照ください。
        前述の関数atan_solverと違い後ろの[5,13]は有りません。)
   
  実行例。

   $ gp atan_solve_check.gp asc_sample.gp
   (起動メッセージ省略)

   ? (asc_sample.gpに記載しているマチンの公式などの定義例が
      実行されて表示されます。)

   ? atan_solver([10,239,515])  (既知の公式を入れてみる)
   8946, 1, +8a(1/10)-a(1/239)-4a(1/515)., [13,101]. (計算結果表示)

   ? atan_check([18,57,239],[12,8,-5])  (既知の公式を入れてみる)
   8933, 1, +12a(1/18)+8a(1/57)-5a(1/239). (計算結果表示)

   ? atan_check([3,4,5],[6,7,8]) (無効な式を入れてみる。)
   NG: 25938, 6.6520634890039804147246471621809574857, +6a(1/3)+7a(1/4)+8a(1/5).

   (終了)
   ? quit
   
*/

/*
関係式のリストの全リストを新形式(2015年以降)から旧形式(2014年以前)に
変換する手順

関係式のリストの全リストは2015年以降データが非常に大きくなっています。
特に6項式はファイルサイズ2Gbyte越えてます。そのため、データ形式を変更し
ています。

   4項式全リスト
     http://www.pluto.ai.kyutech.ac.jp/~matumoto/atan_table_k4.txt.gz
     http://www.pluto.ai.kyutech.ac.jp/~matumoto/atan_table_k4-diff2016.txt.gz
   5項式全リスト
     http://www.pluto.ai.kyutech.ac.jp/~matumoto/atan_table_k5.txt.gz
     http://www.pluto.ai.kyutech.ac.jp/~matumoto/atan_table_k5-diff2016.txt.gz
   6項式全リスト
     http://www.pluto.ai.kyutech.ac.jp/~matumoto/atan_table_k6.txt.gz
     http://www.pluto.ai.kyutech.ac.jp/~matumoto/atan_table_k6-diff2016.txt.gz
   7項式以上 
     http://www.pluto.ai.kyutech.ac.jp/~matumoto/atan_table_k7ov.txt.gz
     http://www.pluto.ai.kyutech.ac.jp/~matumoto/atan_table_k7ov-diff2016.txt.gz
   整数の逆数以外が含まれる公式 
     http://www.pluto.ai.kyutech.ac.jp/~matumoto/atan_table_u-all.txt.gz


新形式 (注意: ] [の間に,(コンマ)が有りません。空白が入っています。)
----
[57,239,682,12943] [44,7,-12,24]
[28,443,1393,11018] [22,2,-5,-10]
[23,182,5118,6072] [17,8,10,5]
----

最初の大[]カッコはa(1/X)内の数字Xを列挙, 二つ目の大カッコ[]はa(1/X)の乗
数です。

旧形式(2014年以前)への変換が必要な場合は本ファイルで定義している関数
atan_solver()/atan_check()を実行してください。但し素数リストの生成には
コストがかかるため、素数リスト有りで生成するプログラムと無しで生成する
プログラムを掲載しています。

旧形式(素数リスト有り)
----
7930, 1, +44a(1/57)+7a(1/239)-12a(1/682)+24a(1/12943)., [5,13,61].
8172, 1, +22a(1/28)+2a(1/443)-5a(1/1393)-10a(1/11018)., [5,157,197].
8554, 1, +17a(1/23)+8a(1/182)+10a(1/5118)+5a(1/6072)., [5,53,373].
----

旧形式(素数リスト無し)
----
7930, 1, +44a(1/57)+7a(1/239)-12a(1/682)+24a(1/12943).
8172, 1, +22a(1/28)+2a(1/443)-5a(1/1393)-10a(1/11018).
8554, 1, +17a(1/23)+8a(1/182)+10a(1/5118)+5a(1/6072).
----

*旧形式(素数リスト有り):

旧形式への変換は以下のとおりです。 コマンドはLinuxを想定。XXをk4,k5,k6,k7ov,u-allに置
き換えてください。(u-allはatan_table_XX-diff20YY.txt.gzは無いです)

中間形式に変換:
   gzip -dc atan_table_XX.txt.gz atan_table_XX-diff20YY.txt.gz | awk '{print "atan_solver("$1");"}' > tmp.gp
      -> tmp.gpを生成

旧形式に変換:
   cat tmp.gp | gp -q atan_solve_check.gp > tmp.txt
      -> 旧形式でtmp.txtに出力。ただしソートされていない。
   sort -n tmp.txt > old-format_XX.txt
      -> 旧形式のold-format_XX.txtというファイルが生成。
   rm tmp.txt tmp.gp
      -> 作業ファイルを削除。  

変換時間:(Intel(R) Xeon(R) CPU E31240 @ 3.30GHz + CentOS7 で測定)

  4項式(k4) -> 1-2数秒
  5項式(k5) -> 3分27秒
  6項式(k6)　-> 5日19時間40分 (79Gのファイルが生成されます)
  7項式以上(k7ov) -> 31分20秒
  逆数以外(u-all)  -> 4時間13分

素数リストで注意点:
巨大な整数の素因数分解の時間短縮のため擬似素数を用いてます。擬似素数を
採用した場合は、素数リストにpseudoと表示されます。

擬似素数が含まれる例。
----
10011, 1, +4a(1/5)-2a(1/478)+a(1/54608383)-a(1/298207603995852)., [pseudo,5,13,41,45697,7273354863109].
----

*旧形式(素数リスト無し)
素数リスト有りで6項式が非常に時間がかかる理由は素数リストの生成に素因数
分解を行ってるからです。素数リストを省略してよければ、以下のコマンドで
中間形式のtmp.gpを生成すれば6項式でも現実的な時間で終わります。

中間形式に変換:
   gzip -dc atan_table_XX.txt.gz  atan_table_XX-diff20YY.txt.gz | awk '{print "atan_check("$1","$2");"}' > tmp.gp
      -> tmp.gpを生成

変換時間:(Intel(R) Xeon(R) CPU E31240 @ 3.30GHz + CentOS7 で測定)
  4項式(k4) ->  1-2秒
  5項式(k5) -> 2分4秒
  6項式(k6)　-> 23時間19分 (46Gのファイルが生成されます)
  7項式以上(k7ov) -> 33分25秒
  逆数以外(u-all)  -> 11時間9分 (遅くなっています)

*PARI/GP実行時間を速くするための補足

gpコマンド実行時に別の端末からtopコマンドを実行してgpのCPU占有率を見て
みてください。上記の処理では常時99%になると思います。90%以下の場合はgp
が作業ファイルディレクトリ/tmpに作るファイルの書き込みに時間がかかって
いる可能性があります。/tmpをRAM-DISKに変更すると改善する可能性がありま
す。CentOS7ではRAM-DISKは/dev/shmです。以下のコマンドを打ち込んで作業ファ
イルディレクトリを変更して top コマンドでCPU占有率を比較してみてくださ
い。

    export GPTMPDIR=/dev/shm 

常時設定する場合は ~/.bashrcに加筆すれば良いです。

以上

*/


/**************************************************************************

                    以下プログラム本体

**************************************************************************/


/* 
 atan_solver制御パラメーター

 0の場合、うまく立式できなかった場合でも出力する。 
 2の場合、うまく立式できた場合のみ出力する。
*/
flag_zero=0;

/*
  公式の評価値 n の表示形式
  0: 円周率10000桁を計算するのに必要なテイラー級数の総和。
  1  1/log(X)の総和。
*/
flag_logsum=0;

/* 2015/1/18, atan_solver_aux()より分離 */
print_logsum(X, K) =
{
    \\改行無し
    local(i, ret, xlen, num);
    \\print("print_logsum: X= ", X);
    xlen = length(X);
    num = 0.0;
    for(i=1,xlen,
	    if(K[i]==0, next);
	    num += 1 / log(X[i]);
	);

    if (flag_logsum==0, 
	    print1(round(log(10)*5000*num)),
	    printf("%.8f", num)
    );
}

/*
   係数行列の符号計算の関数。詳細は後述の説明参照
*/
atan_solver(X) =
{
    local(ret);
    ret=-1;
    ret = atan_solver0(X, []); \\ ステルマーのオリジナル+整数の逆数以外に拡張。
    \\ret = atan_solver1(X, []); \\ pari/gpのfactor()を使う。
    return(ret);
}


/* 2015/1/18追加
   この関数はatan_solver_aux()が使っている。

flag_print_plist

      1: 素数リストを表示する
      0: 素数リストを表示しない
*/
flag_print_plist=1;

print_plist(P) = 
{
    local(i, plen, s_p);
    if(flag_print_plist==0, return());
    plen = length(P);
    s_p =vecsort(P);
    print1("., [");
    for(i =1,  plen,
	if(s_p[i] == -1, print1("pseudo"), print1(s_p[i]));
	    if(i < plen, print1(",")););
    print1("]");
}


/***************************************************

      以下の変数、 prime_collect()の制御に使う変数。
      この関数はatan_solver0()が使っている。

****************************************************/

/*
  素数収集の時に極端に大きいX[i]に対して
   0: 常にpari/gpのfactor()を使う。
   1: 手抜き素因数分解を行う。擬似素数が含まれる可能性が有る。

注意事項: 0 にすると手抜き素因数分解を行いませんが、非常に遅くなります。
2015年1月以降の6項式全リスト(atan_table_k6.txt.gz)に対して0で処理すると
数ヶ月かかると思います。

制限事項: 手抜き素因数分解は、整数の逆数以外の関係式の場合、
早々にあきらめて擬似素数と評価する可能性有り。atan(y/x)の
第一項の(y^2+x^2)がある程度小さい値である事を想定してプログラム
しているが、整数の逆数以外の関係式の場合、第一項から大きい値の
場合がある。
*/
flag_tenuki_factor=1;

/*
  手抜き素因数分解の閾値。flag_tenuki_factor=1の場合、以下の数値
　より大きい時は擬似素数が含まれる可能性が有る方法を取る。

  atan(y/x)で、 (x^2 + y^2) > Normal_factor_limit  のときに手抜き。
*/
Normal_factor_limit = 10^12;
\\Normal_factor_limit = 10^4;
\\Normal_factor_limit = 100;

/*
 擬似素数と評価する時は、factor()の第二引数で素数の範囲を制限して評価する。
 制限された素数の範囲内で割れなければ擬似素数の評価する。
 Ref: http://pari.math.u-bordeaux.fr/dochtml/html.stable/Arithmetic_functions.html

 第二引数で渡す素数の範囲

 通常gpを起動した時に表示されるprimelimitの値を書いておけばよいと思います。
*/
Factor_Arg2_lim = 500000; 

\\default(primelimit, Factor_Arg2_lim);
\\if(Factor_Arg2_lim > default(primelimit), Factor_Arg2_lim = default(primelimit));
\\print(default(primelimit));

/*
  手抜き素因数分解flag_tenuki_factor=1を使う時、

   0: 正常に素因数分解を行った。手抜き素因数分解を実施しなかった。
   1: 手抜き素因数分解が行われた。擬似素数が含まれる可能性がある。

　この変数はprime_collect()が呼ばれるたびに0にリセットされ、
  手抜き素因数分解が実施されると1になる。
*/
flag_have_pseudo_prime=0;


/*********************************************************************
            
                    関数定義部
 
*********************************************************************/

/*
  atanのk概算
  2015/1/18, atan_solver_aux()より分離
*/

calc_atan_k(X, K) =
{
    local(s,k,i,xlen);
    xlen = length(X);
    s = 0;
    for(i=1,xlen, s+=atan(1/X[i])*K[i]);
    \\ Pi/4の何倍か計算する。
    k = s / (Pi/4); 
    return(k);
}

/*
  atanの公式表示
  2015/1/18, atan_solver_aux()より分離
*/
print_alist(X, K) = 
{
    local(i, xlen);
    xlen = length(X);
    for(i=1,xlen,
	if(K[i]>=0,print1("+"));
	if(K[i]<0,print1("-"));
	if(abs(K[i])==1, print1("a(",1/X[i],")"),
	               print1(abs(K[i]),"a(",1/X[i],")"));
    );
}

/*
  素数収集
*/
flag_collect_prime_debug=0;
flag_collect_prime_tenuki_debug=0;
collect_prime(X) = 
{
    local(i, j, xlen, P, plen, p, x2n1, Y);
    flag_have_pseudo_prime=0;
    P = [];
    xlen = length(X);
    if(xlen == 0, return (P));
    Y = [];
    for(i = 1, xlen, 
	x2n1= numerator(X[i])^2 + denominator(X[i])^2;\\整数の逆数以外に拡張(2015/1/5)
	if(x2n1%2==0, x2n1=x2n1/2);\\素数2で割る。

        \\既知の素数で割る。
        plen = length(P);
        for(j = 1, plen, while(x2n1%P[j]==0, x2n1=x2n1/P[j]));
        if(flag_tenuki_factor==1 && x2n1 > Normal_factor_limit, Y = concat(Y, x2n1); next);
	p = factor(x2n1);
	p = p[,1]~;
	P = setunion(P, p);
    );

    if(length(Y)>0, 
	    if(flag_collect_prime_debug==1,
		    print("# in collect_prime()");
		    print("# using tenuki, X=",X);
		    print("#               Y=",Y);
		    print("#               P=",P);
	    );
	    return(collect_prime_tenuki(Y, P))
    );
    return(P);  
}

/*
  素数収集, 手抜き版
*/

collect_prime_tenuki(X,P) = 
{
    local(i, j, g, xlen, p, y1, Y, ylen, Z, zlen, flag, plen);

    Y=[];    
    xlen = length(X);
    plen = length(P);
    for(i = 1, xlen, Y=concat(Y,X[i]); for(j = 1, plen, while(Y[i]%P[j]==0, Y[i]=Y[i]/P[j])));

    Z=[];
    for(i = 1, xlen, if(Y[i]>1, Z=setunion(Z,[Y[i]])));
    Y=Z;

    ylen = length(Y);
    if(ylen == 0, return (P));

    y1 = Y[1];
    Z=[];
    for(i = 2, ylen, Z = concat(Z, Y[i]));

    \\ 制限値以下は通常の素因数分解
    if(y1 <= Normal_factor_limit,
	p = factor(y1);
	p = p[,1]~;
	P = setunion(P, p);
	return (collect_prime_tenuki(Z,P));
    );

    \\ GCDなどする素因数分解
    Z=[];
    flag = 0;
    for(i=2, ylen,
	    g=gcd(Y[i],y1);
	    if(g == 1, Z=setunion(Z, [Y[i]]); next);
	    flag++;
            Z = setunion(Z, [g]);
            Z = setunion(Z, [Y[i]/g]);
            Z = setunion(Z, [y1/g]);
    );

    if(flag > 0, return (collect_prime_tenuki(Z,P)));

    \\ すなおに素因数分解すると時間かかるので制限付き

    if(flag_collect_prime_tenuki_debug==1,
	    print1("# In collect_prime_tenuki(): factor of a very large number : ",y1));
    p = factor(y1, Factor_Arg2_lim);
    p = p[,1]~;

    \\ 因数分解失敗。擬似素数のフラグを立てる。
    if(p[1] == y1, 
	    flag_have_pseudo_prime=1;
	    if(flag_collect_prime_tenuki_debug==1,
		    print(" : fail, using pseudo prime."));
	    P = setunion(P, [y1]);
	    return(collect_prime_tenuki(Z,P))
    );

    \\ [制限付き]の素因数分解成功。

    plen = length(p);
    flag = 0;
    for(i=1, plen, 
	if(p[i] < Factor_Arg2_lim,
		P = setunion(P, [p[i]]), \\ 制限値以下は素数として登録
		Z = setunion(Z, [p[i]]);flag=1 \\制限値より大きいのは再度試す。
        )
    );
    if(flag_collect_prime_tenuki_debug==1,
	    if(flag==0, print(" : ok."), print(" : retry.")));
    return(collect_prime_tenuki(Z,P));
}

/* 
  atan_solver() 関係の関数群

  Ref: http://www.jjj.de/arctan/arctanpage.html
  Ref: http://cm.hit-u.ac.jp/~kobayashi/topics/arctan.pdf  
*/

atan_solver_aux(X, P, M) =
{
    local(K, K1, klen, xlen, stat, k, i, j, flag_have_zero, ret_stat);
    ret_stat=0;
    stat = 0;
    \\ 係数を算出。
    K=matkerint(M);
    klen = length(K);
    xlen = length(X);
    \\ Kが空行列だった場合、ゼロ行列を作る。
    if(klen == 0,
	    ret_stat=1;
	    if(flag_zero > 0, return(ret_stat));
   	    K=matrix(xlen,1);
            stat = 1
    );
    klen= matsize(K);\\二重ベクトルの場合がある
    if(klen[2] > 1, stat=2; if(flag_zero > 0,  return(1)));
    
    for(i=1, klen[2],
	K1 = K[,i]~;
	flag_have_zero=0; \\ Kにゼロが入っているか否か
	if(setsearch(vecsort(K1), 0) > 0, flag_have_zero=1; ret_stat=1);
	

        \\ 組み込み関数atan()を使った調査
	k = calc_atan_k(X, K1);
	if(k < 0, k = -k; K1=-K1);\\ atanの式をプラスに修正する。
	\\表示
	if(abs(k-round(k)) > 0.01 , print1("NG: "), k=round(k));

	if(flag_zero >= 2 && k==0, next);

	\\ゼロ解(0 * PI/4=)の場合は、先頭の符号を+にする。
	if(k==0, if(K1[1]<0, K1=-K1));

	if(stat == 1, print1("Z: "));\\k=0
	if(stat == 2, print1("V: "));\\2重ベクトル
	if(stat == 0 && flag_have_zero == 1, print1("z: "));\\k>0だけどatan(1/x)にゼロ倍が含まれる

        print_logsum(X, K1);
	print1(", ",k,", ");
	print_alist(X, K1);
	print_plist(P);
	print(".");
     );
     return(ret_stat);
}

/* 
- atan_solver0(X,P) 

 連立方程式の符号確定に、ステルマーのオリジナルの方法を使う。
 Pは空でも動くが素因数分解するので効率が悪い。引数で与えるのがベター。
 但し、Xに分数が可能なように拡張。Xは約分済みであること。

 一部自作アルゴリズムを用いているので、結果が誤っていることを想定してください。
*/

/* 
  atan_solver0_sign(X,P)
  2015/1/18, atan_solver_aux()より符号確定部分離
 */
atan_solver0_sign(X,P) = 
{
    local(xlen, plen, M, i, j, x2n1, flag, flag_logical_check, tmp1, tmp2);

    xlen = length(X);
    plen = length(P);
    M = matrix(plen, xlen); \\係数の行列。
    
    \\ 行列の作成
    for(i=1,xlen,
        \\ numerator()   : 分子, 但しX[]は逆数のため、実際は分母
        \\ denominator() : 分母, 但しX[]は逆数のため、実際は分子
	x2n1=denominator(X[i])^2 + numerator(X[i])^2;
	for(j=1,plen,
	    while(x2n1%P[j]==0,
		M[j,i]=M[j,i]+1;
		x2n1=x2n1/P[j]
	    );
	);
	if(x2n1 > 2, print("ERROR: atan_solover0: Logical Error: X = ", X); return([]));
    );
    for(i=1, plen,
	flag=0;
	for(j=1, xlen,
	    if(M[i,j] == 0, next);
	    if(flag==0, flag=j; next);
	    if(flag > 0, 
		    flag_logical_check=0; \\ アルゴリズムが自作のため、反例が無いか確認
		    \\ numerator()   : 分子, 但しX[]は逆数のため、実際は分母
                    \\ denominator() : 分母, 但しX[]は逆数のため、実際は分子
		    tmp1 = numerator(X[j])*denominator(X[flag]);
		    tmp2 = numerator(X[flag])*denominator(X[j]);
		    if((tmp1+tmp2)%P[i] == 0, flag_logical_check++; M[i,j]=-M[i,j]);
		    if((tmp1-tmp2)%P[i] == 0, flag_logical_check++);
	            if(flag_logical_check != 1, print("ERROR: atan_solover0: Logical Error: X = ", X); return([]))
 	    )
       )
    );
    return (M);
}

atan_solver0(X,P) = 
{
    local(M);

    if(length(X) == 0, print("ERROR: atan_solver0: Bad Input X=",X, ", P=", P); return(0));
    \\ 素数表。空でも自動的につくる。
    if(length(P) == 0, P = collect_prime(X));

    M = atan_solver0_sign(X, P);

    if(length(M)==0, return(0));

    if(flag_have_pseudo_prime==1, P=concat(-1, P));\\擬似素数が有る場合は素数リストに[-1]を追加。
    return (atan_solver_aux(X, P, M));
}

/*
- atan_solver1(X, P) 

連立方程式の符号確定にpari/gpのfactor()をもちいてガウス整数の素因数分解
を行なう。Pは空でも、与えても効率には影響しない。常に素因数分解する。X
は分数でも可。但し約分済みであること。

自作のアルゴリズムが混ざってないため、念のための実行時はこれを用いる。

ガウス整数の因数分解を毎回行っているが、高速化を目指すなら、事前に因数
分解を行い、表にし方がよいが、面倒なので省略。

実際に用いているC++言語版のプログラムでは表にしている。
*/
atan_solver1(X, P) =
{
    local(xlen, plen, M, i, j, pp, n, k, flag, b);

    xlen = length(X);
    plen = length(P);\\ 素数表。空でも自動的につくる。
    M = matrix(xlen,plen); \\係数の行列。プログラムの都合でMは転置行列を作っている。
    for(i=1,xlen,
	\\ガウス整数の因数分解
	b = factor(numerator(X[i])+denominator(X[i])*I);
	\\行列M作成
	for(j=1,length(b~),
            pp = b[j,1];
	    n = norm(pp);
	    if(n==1 || n==2, next); \\ Gauss整数の単数と1+iを無視
	    plen = length(P);
	    flag=0;
	    for(k=1, plen, if(P[k] == n, flag=k; break));
            if(flag==0, \\素数表Pに無い素数の場合P,Mを拡張する。
		flag = plen + 1;
		P = concat(P, [n]);
		M = concat(M, matrix(xlen,1));
	    );
            \\符号判定
	    if(real(pp) < imag(pp), b[j,1] *= -I);
	    M[i,flag]=b[j,2] * sign(imag(b[j,1]));
    ););
    M = M~;
    return (atan_solver_aux(X, P, M));
}

/*
- atan_solver2(X, P, S)

連立方程式の符号を引数Sで与える。Pは必ず与えること。この関数のみ引数の
数が違う。

gp2cでC言語に変換して使うためのインターフェース用。

検索用の本体プログラムは非公開です。面倒な素因数分解等は、C言語側で表を
用いて行っている。

Sは本来行列であるがC言語のインターフェースで用いているため、引数を渡すのが
面倒なのでベクトルにしている。
*/
atan_solver2(X, P, S) = 
{
    local(xlen, plen, slen, M, i, j);

    xlen = length(X);
    plen = length(P);
    slen = length(S);
    if(xlen * plen != slen || xlen ==0 || plen == 0,
	print("ERROR: atan_solver2: Bad Input X=",X, ", P=", P, ", S=", S); return(0));

    M = matrix(plen, xlen); \\係数の行列。

    \\ Sは本来行列だが、面倒なので、ベクトルで渡して、行列に読み替え。
    for(i=1, length(P),
	for(j=1, xlen,
	    M[i,j] = S[(j-1)*plen + i]
	)
    );
    return(atan_solver_aux(X, P, M));
}

/*
既知の公式のチェックをするプログラム

 Ref: http://mathworld.wolfram.com/Machin-LikeFormulas.html
 http://www.pluto.ai.kyutech.ac.jp/~matumoto/dvi/pi.pdf, page 33.

 2重ベクトルでも正常と判断されるので注意。
*/

/*
 2015/1/18追加。 (実質的な)通分関数
 絶対値は変わっても問題なく、偏角のみ情報が必要なため。
*/
my_gcd(mul_g) = 
{
    local(g);
    \\return(mul_g);\\ debug;

    g = gcd(abs(real(mul_g)), abs(imag(mul_g)));
    if(g > 1, return(mul_g / g), return(mul_g));
}

/*
  atan_check_aux_mul()
  複素数を使った調査。 atan_check(X, K) より分離
*/
atan_check_aux_mul(X, K, k) = 
{
    local(mul, xlen, i, mul0, tmp, flag_ok);
    xlen = length(X);
    mul = 1;
    for(i=1, xlen,
	tmp = numerator(X[i])-denominator(X[i])*I;
	mul0=tmp^abs(K[i]);
	/*
	j = abs(K[i]);
	mul0=1;
	while(j>0, 
	    if(j%2==1, mul0 = my_gcd(mul0 * tmp);j=j-1);
		j=j/2;
		if(j>0, tmp = my_gcd(tmp^2))
	);*/
	mul0 = my_gcd(mul0);
	\\ mul = mul * conj(mul0)は, mul = mul / mul0 と同じ
	if(K[i] > 0, mul = mul*mul0, mul = mul*conj(mul0));
            mul = my_gcd(mul);			   
    );
    \\print("mul_log=",log(abs(mul))/log(10));
    flag_ok = 0;
    if(k % 4 == 1 && real(mul) == -imag(mul), flag_ok=1);
    if(k % 4 == 2 && real(mul) == 0, flag_ok=1);
    if(k % 4 == 3 && real(mul) == imag(mul), flag_ok=1);
    if(k % 4 == 0 && imag(mul) == 0, flag_ok=1);
    if(flag_ok == 0, return(1));
    return(0);
}

/*
 2015/1/18追加。
 ばかでっかい式は複素数で計算すると巨大になりメモリ不足を起こすので閾値以上は
 atan_solver0()のエンジンを呼び出して演算行列を作り
 それを連立方程式として解くのではなく、Kを代入して0になる事を確認。
*/
atan_check_aux_stoermers_sign(X, K) = 
{
    local(M, P, xlen, plen, i, j, s);
    P = collect_prime(X);
    M = atan_solver0_sign(X, P);    
    if(length(M)==0, return(1));
    xlen = length(X);
    plen = length(P);

    if(flag_logsum2_debug==1, print("M=",M));

    for(i=1, plen, for(j=1, xlen, M[i, j] *= K[j]));
    
    for(i=1, plen, s=0; for(j=1, xlen, s+= M[i, j]); if(s!=0, return(1)));
    return(0);
}

\\ ばかでっかい式の評価変更のdebug
flag_logsum2_debug =0; \\ disable
\\flag_logsum2_debug =1; \\ enable

atan_check(X, K) =  
{
    local(xlen, klen, logsum2, flag_have_zero, k, i, ret_stat, g, tmp, flag_bignum);

    ret_stat = 0;
    xlen = length(X);
    klen = length(K);
    \\K[klen]=1; \\ noize

    if(xlen != klen, 
	    print("ERROR: Bad Input X=",X, ", K=", K); return (1)
    );

    \\ 最大公約数で割る。多くの場合不要ですが。
    g = abs(gcd(K));
    if(g > 1,K = K/g);

    \\ 事前にパラメータ計算
    flag_have_zero=0; \\ Kにゼロが入っているか否か
    if(setsearch(vecsort(K), 0) > 0, flag_have_zero=1);
          
    \\ 組み込み関数atan()を使った調査
    k = calc_atan_k(X, K);
    if(k < 0, k = -k; K = -K);\\ atanの式をプラスに修正する。
    if(abs(k-round(k)) > 0.01, print1("NG: "); ret_stat=1, k=round(k));

    \\この計算式で大体atan_check_aux_mul()の内部計算で使う複素数の絶対値になります。
    \\なお、実際はmy_gcd()でしつこく通分しているので、もうちょっと小さい値になってます。
    logsum2 = 0; \\ あまりにでかい式はあきらめる。
    for(i=1,xlen,logsum2 += log(X[i]) * abs(K[i]));
    logsum2 /= log(10);

   if(flag_logsum2_debug==1,
	printf("%.3f, ", logsum2);
	print(X,K));

   \\ 複素数を使った調査
   flag_bignum=0;
   if(ret_stat == 0, 
       \\ if(100000 < logsum2, \\ あまりにでかい式はあきらめて、別の評価関数を使う
       if(10000 < logsum2, \\ あまりにでかい式はあきらめて、別の評価関数を使う
		if(flag_logsum2_debug==1, 
                     print("big num, call atan_solver(): num = ",logsum2));
 		     flag_bignum=1;
		     ret_stat = atan_check_aux_stoermers_sign(X, K),
		\\ else
		     ret_stat = atan_check_aux_mul(X, K, k)
	);
	if(ret_stat == 1, print1("NG: "));
   );

   \\ k>0だけどatan(1/x)にゼロ倍が含まれる
   if(flag_have_zero == 1 && ret_stat == 0, 
        print1("z: "); ret_stat=1);

   print_logsum(X, K); 
   print1(", ",k,", ");
   print_alist(X, K);
   print(".");
   \\if(flag_bignum==1, print("# big"));
   return (ret_stat);
}

/* EOF */
