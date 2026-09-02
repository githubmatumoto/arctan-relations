\\ -*- sh -*-
/*
  ライセンス: GPL
  作者: 松元隆二 (matsumoto(AT)tech-i.kyutech.ac.jp)
  最終更新日: 2014/12/13

  1. arctan関係式の係数を計算するプログラム
       - atan_solver()

  2. 既知の公式のチェックをするプログラム
       - atan_check()
        
  PARI/GPが必要です。
    http://pari.math.u-bordeaux.fr/

  使い方はファイルの最後を見てください。
*/

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

/* 
  atan_solver() 関係の関数群

  Ref: http://www.jjj.de/arctan/arctanpage.html
  Ref: http://cm.hit-u.ac.jp/~kobayashi/topics/arctan.pdf  
*/

atan_solver_aux(X, P, M) =
{
    local(K, klen, xlen, stat, s, k, i, j, logatan, flag_have_zero, ret_stat, s_p, s_p_l);
    ret_stat = 0;
    stat = 0;
    \\ 係数を算出。
    K=matkerint(M);
    klen = length(K);
    xlen = length(X);
    \\ Kが空行列だった場合、ゼロ行列を作る。
    if(klen == 0,
	    if(flag_zero > 0, return(ret_stat));K=matrix(xlen,1); stat = 1,
	    ret_stat = 1;
	    
    );
    klen= matsize(K);\\二重ベクトルの場合がある
    if(klen[2] > 1, stat=2; if(flag_zero > 0,  return(ret_stat)));
    
    for(i=1, klen[2],
	logatan=0.0; \\級数の収束の速さ
	\\ atanの式をプラスに修正する。
	s=0;
	flag_have_zero=0;
	for(j=1,xlen,
	    if(K[j,i]==0, flag_have_zero=1; next);
	    s+=atan(1/X[j])*K[j,i];
	    logatan += 1 / log(X[j]);
	);

	if(s<0, s = -s; K=-K);
	\\ Pi/4の何倍か計算する。
	k = s / (Pi/4); 
	\\表示
	if(abs(k-round(k)) > 0.01 , print1("NG: "), k=round(k));

	if(flag_zero >= 2 && k==0, next);

	\\ゼロ解の場合は、先頭の符号を+にする。
	if(k==0, if(K[1,i]<0, for(j=1,xlen,K[j,i]=-K[j,i])));

	if(stat == 1, print1("Z: "));\\k=0
	if(stat == 2, print1("V: "));\\2重ベクトル
	if(stat == 0 && flag_have_zero == 1, print1("z: "));\\k>0だけどatan(1/x)にゼロ倍が含まれる

	if (flag_logsum==0, 
		print1(round(log(10)*5000*logatan)),
		printf("%.8f", logatan,k));
	print1(", ",k,", ");

	for(j=1,xlen,
	    if(K[j,i]>=0,print1("+"));
	    if(K[j,i]<0,print1("-"));
	    if(abs(K[j,i])==1, print1("a(",1/X[j],")"),
		               print1(abs(K[j,i]),"a(",1/X[j],")"));
	);

	s_p =vecsort(P);
	s_p_l = length(s_p);
	print1("., [");
	for(j =1,  s_p_l,
	    print1(s_p[j]);
	    if(j < s_p_l, print1(",")););
	    
	print("].");
	
     );
     return(ret_stat);
}

/*
  係数計算のトップ関数。詳細は後述の説明参照
*/
atan_solver(X) =
{
    return (atan_solver1(X, []));
}

/*
   説明省略。

  atan_solver()では省略していた、共通の素数を直接指定する場合。
  なお、共通の素数を指定しても早くなるわけではない。。

　ガウス整数の因数分解を毎回行っているが、高速化を目指すなら、
　事前に因数分解を行い、表にし方がよいが、面倒なので省略。

　実際に用いているC++言語版のプログラムでは表にしている。
*/
atan_solver1(X, P) =  \\ Pは空でも動く。
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

\\ 既知の公式のチェックをするプログラム
\\
\\ Ref: http://mathworld.wolfram.com/Machin-LikeFormulas.html
\\ http://www.pluto.ai.kyutech.ac.jp/~matumoto/dvi/pi.pdf, page 33.
\\
\\ 二重ベクトルでも正常と判断されるので注意。
\\
atan_check(X, K) =  
{
    local(xlen, klen, logatan, flag_have_zero, s, g,i, mul, mul0, tmp, flag_ok, ret_stat);

    ret_stat = 0;
    xlen = length(X);
    klen = length(K);

    if(xlen != klen, 
	    print("ERROR: Bad Input X=",X, ", K=", K); return (1)
    );

    g = abs(gcd(K));
    if(g > 1,K = K/g);

    logatan = 0.0;\\級数の収束の速さ
    s = 0; \\ atanの式をプラスに修正する。
    flag_have_zero=0;
    for(i=1,xlen,
	    if(K[i]==0, flag_have_zero=1; next);
	    s+=atan(1/X[i])*K[i];
	    logatan += 1 / log(X[i]);
	);
   if(s<0, s = -s; K=-K);
   k = s / (Pi/4); 
   \\表示
   if(abs(k-round(k)) > 0.01, print1("NG: "); ret_stat=1, k=round(k));

   if(ret_stat == 0,
	   mul = 1;
	   for(i=1, xlen,
	       mul0=1;
	       tmp = numerator(X[i])-denominator(X[i])*I;
	       for(j=1, abs(K[i]),mul0 = mul0 * tmp);
	       if(K[i] > 0, mul = mul*mul0, mul = mul/mul0);
	   );
           \\print("mul=",mul);
	   flag_ok = 0;
	   if(k % 4 == 1 && real(mul) == -imag(mul), flag_ok=1);
	   if(k % 4 == 2 && real(mul) == 0, flag_ok=1);
	   if(k % 4 == 3 && real(mul) == imag(mul), flag_ok=1);
	   if(k % 4 == 0 && imag(mul) == 0, flag_ok=1);
   
           if(flag_ok == 0, print1("NG: ");ret_stat=1);
    );
   if(flag_have_zero == 1 && ret_stat == 0, print1("z: "); ret_stat=1);\\k>0だけどatan(1/x)にゼロ倍が含まれる
   if (flag_logsum==0, 
	print1(round(log(10)*5000*logatan)),
	printf("%.8f", logatan,k));

   print1(", ",k,", ");

   for(i=1,xlen,
        if(K[i]>=0,print1("+"));
        if(K[i]<0,print1("-"));
        if(abs(K[i])==1, print1("a(",1/X[i],")"),
               print1(abs(K[i]),"a(",1/X[i],")"));
    );
    print(".");
    return (ret_stat);
}

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
        (出力フォーマットは atan_table.txt参照ください。)

  とすれば、算出されます。

  2. 既知の公式のチェックをするプログラム

　関数atan_check()で既知の公式のチェックをします。ただし、二重ベクトルの式も正常と認識します。

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
       (出力フォーマットは atan_table.txt参照ください。
        前述の関数atan_solverと違い後ろの[5,13]は有りません。)
   
  実行例。

   $ gp atan_solve_check.gp
   (起動メッセージ省略)

   ? (このファイルの最後に記載しているマチンの公式などの定義例が実行されて表示されます。)

   ? atan_solver([10,239,515])  (既知の公式を入れてみる)
   8946, 1, +8a(1/10)-a(1/239)-4a(1/515)., [13,101]. (計算結果表示)

   ? atan_check([18,57,239],[12,8,-5])  (既知の公式を入れてみる)
   8933, 1, +12a(1/18)+8a(1/57)-5a(1/239). (計算結果表示)

   ? atan_check([3,4,5],[6,7,8]) (無効な式を入れてみる。)
   NG: 25938, 6.6520634890039804147246471621809574857, +6a(1/3)+7a(1/4)+8a(1/5).

   (終了)
   ? quit
   
*/

\\
print("# atan_solver")

\\マチンの公式
atan_solver([5,239]); 

\\ガウスの公式
atan_solver([18,57,239]); 

\\ガウスの公式
atan_solver([18,57,239]); 

\\ オイラーの公式。整数の逆数で無い式。
atan_solver([7,79/3]); 

\\ 無効な式。先頭にZが表示されます。
atan_solver([8,18]); 

\\ 逆数で無い式で無効な式。
atan_solver([8,18/5]); 

\\2重ベクトルの式。2式以上解が表示され先頭にVが表示されます。
atan_solver([2,3,7]);

\\係数kにゼロが含まれる式。先頭にzが表示されます。
atan_solver([5,122,229,239,1252,246303]);

\\K=0 (pi/4*K = .. )の式。式は有効ですが、円周率の計算には使えません。
atan_solver([4,8,268]);

\\
print("# atan_check")

\\マチンの公式
atan_check([5,239],[4, -1]); 

\\ オイラーの公式。整数の逆数で無い式。
atan_check([7,79/3],[5,2]); 

\\ 2重ベクトルの式。正常と判断されます。
atan_check([2,3,7],[1,3,1]);

\\K=0 (pi/4*K = .. )の式。式は有効ですが、円周率の計算には使えません。
atan_check([4,8,268],[1,-2,1]);

\\無効な式。NGと表示されます。
atan_check([10,20],[1,2]);

\\係数kにゼロが含まれる式。先頭にzが表示されます。
atan_check([5,122,229,239,1252,246303],[4,0,0,-1,0,0]);

print ("");

/* EOF */
