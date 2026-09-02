\\ -*- gp-script -*-
/*
  ライセンス: GPL
  作者: 松元隆二 (matsumoto(AT)tech-i.kyutech.ac.jp)
  last updated : 2022/06/08 (JST)

  PARI/GPが必要です。
    http://pari.math.u-bordeaux.fr/

2022/6/8: 公開版作成

*/

/*

Arctan関係式の探索に用いる「x^2+1」の素因数分解の表を作成するプログラム
です。このプログラムは簡易的なものです。本格的に作る場合はC言語などで実
行した方がよいです。

補足: 他の自作プログラムの関係でnは逆数(1/n)として扱います.

  1. n^2+1を素因数分解する。 
     - bun(n):  

  2.  分母(n)^2+分子(n)^2を計算する。
     - bun_u(n):  

  3. 上記bun_u(n)に追加で最大の素因数引数を指定。
     素因数がlimを越えたら何も表示しない。

     - bun_u_lim(n, lim)

  上記いずれも出力は

     n:素因数分解の結果を「,」で並べて最後は「.」

 になります。

例1:
     入力
     bun(1000);

     計算
     1000^2+1 = 101 * 9901

     出力
     1000:101,9901.

例2:
     入力
     bun_u(100/33);

     計算
     100^2+33^2 = 13 * 853

     出力
     100/33:13,853.

例3: 最大の素因数を指定。

     入力
     bun_u_lim(100/33,100);

     計算
     100^2+33^2 = 13 * 853

     出力
     (素因数が100を越えたので何も表示しない。)

例4: pari/gpが通分するので、既約でなくてもいい。

     入力
     bun_u(55/5);

     出力
     11:2,61.

例5: とりあえず素因数分解の表を作りたかったら以下のようにすればよいです。
(Linuxなどのコマンドラインが使えるOSで実行ください)

Nで表の作成範囲を指定してください。

   echo "N=10^4;for(i=2,N,bun(i));"  | gp -q ./term-dat.gp

最大の素数を指定したい場合はPで最大の素数を指定してください。

   echo "N=10^4;P=65536;for(i=2,N,bun_u_lim(i,P));"  | gp -q ./term-dat.gp

分数も含める場合は、以下のようにしてください。

   echo "N=10^3;P=140;for(i=2,N,for(j=1,i-1,if(gcd(i,j)==1,bun_u_lim(i/j,P))));" | gp -q term-dat.gp

以上です。
*/

\\ 整数の逆数のみ。
bun(n) = 
{
	local(v,i,j,s);	
	v = factor(n^2+1);\\因数分解
	s = matsize(v);\\大きさを取得
	print1(n);
	i=1;
	while(i<=s[1],
		j=1;
		while(j<=v[i,2],
			if(i==1 & j==1, print1(":"),  print1(","));
			print1(v[i,1]);
			j++;
			);
		i++;
	);
	print(".");
}

\\ 整数の逆数以外でも可

bun_u(n) = 
{
	local(v,i,j,s);	

	\\ numerator()   : 分子
	\\ denominator() : 分母
	v = factor(numerator(n)^2+denominator(n)^2);\\因数分解

	s = matsize(v);\\大きさを取得
	print1(n);
	i=1;
	while(i<=s[1],
		j=1;
		while(j<=v[i,2],
			if(i==1 & j==1, print1(":"),  print1(","));
			print1(v[i,1]);
			j++;
			);
		i++;
	);
	print(".");
}

\\ 素因数がlimより大きい場合は何も表示しない。
bun_u_lim(n,lim) = 
{
	local(v,i,j,s);	

	\\ numerator()   : 分子
	\\ denominator() : 分母
	v = factor(numerator(n)^2+denominator(n)^2);\\因数分解

	s = matsize(v);\\大きさを取得

	\\ limitより大きいなら何も表示しない。
	if(v[s[1],1] >lim, return (0));
	
	print1(n);
	i=1;
	while(i<=s[1],
		j=1;
		while(j<=v[i,2],
			if(i==1 & j==1, print1(":"),  print1(","));
			print1(v[i,1]);
			j++;
			);
		i++;
	);
	print(".");
	return (1);
}

\\EOF
