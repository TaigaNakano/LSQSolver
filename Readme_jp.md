# LSQSolver

 [<u>English version</u>](https://github.com/TaigaNakano/LSQSolver/blob/main/Readme.md)

.NET向け密行列最小二乗ソルバーです。

以下の問題をサポートします。

- 過決定問題
- 劣決定問題
- ランク落ち問題
- 最小ノルム解

SVDを使用せず、列ピボット付きQR分解（CPQR）に基づいて解を求めます。

GitHub リポジトリ: <https://github.com/TaigaNakano/LSQSolver>

---

## 特徴

- 列ピボット付きQR分解（CPQR）
- 自動ランク判定
- 最小ノルム解
- ランク落ち問題への対応
- Parallel.Forによる並列化
- overwriteモードによる省メモリ動作
- QR中間結果の保存機能
- 外部ライブラリ不要

---

## 開発の背景

多くの数値計算ライブラリは、SVDや擬似逆行列を用いて最小二乗問題を解いています。

LSQSolverは、実用上十分な精度を維持しながら、CPQRとCholesky分解を利用してランク落ち問題や劣決定問題に対応することを目的として開発されています。

---

## インストール

```bash
dotnet add package LSQSolver
```

## 使用方法

### モジュールのインポート

```csharp
using LSQSolver;
using static LSQSolver.LSQSolver;
```

### 行列の作成

```csharp
double[][] data =
{
    new double[] { 1, 2 },
    new double[] { 3, 4 },
    new double[] { 5, 6 }
};

var A = new MatrixObject(data);
```

### 最小二乗問題を解く

```csharp
double[] b = { 7, 8, 9 };
var result = Solve(A, b); // or LSQSolver.LSQSolver.Solve(A, b);
```

### 解法 (Solve Method)

```csharp
var result = Solve(
    A,
    b,
    overwrite: true,
    store_intermediates: false,
    rank_tolerance: 2.22044604925032e-16,
    check_finite: true);
```

### パラメータ

- `A`: 係数行列。
- `b`: 右辺ベクトル。
- `overwrite`: `true` の場合、`A` と `b` を上書きしてメモリ割当を削減します。
- `store_intermediates`: `true` の場合、`R`、`Qᵀb`、およびピボット情報を保存します。
- `rank_tolerance`: 数値ランク判定に使う相対許容値。
- `check_finite`: `true` の場合、`A` と `b` に `NaN` や `Infinity` が含まれていないかを検査します。

### 戻り値：　LSQSolverResult

`Solve()` は `LSQSolverResult` オブジェクトを返します。

`LSQSolverResult` には、計算された解、ソルバーのステータス、基本的な診断情報、および必要に応じて保存された QR 分解の中間情報が含まれます。

`Solution` を使用する前に、必ず `Status` を確認してください。

| プロパティ | 型 | 説明 |
| --- | --- | --- |
| `Status` | `LSQSolverStatus` | 解法の実行結果を表すステータスコードです。`Solution` を使用する前に確認してください。 |
| `Solution` | `double[]` | 計算された最小二乗解ベクトルです。ランク落ち問題や劣決定問題では、解法が成功した場合、最小ノルム解が返されます。 |
| `Rows` | `int` | 入力行列 `A` の行数です。 |
| `Cols` | `int` | 入力行列 `A` の列数です。 |
| `Rank` | `int` | 列ピボット付き QR 分解の過程で推定された数値ランクです。 |
| `ResidualNorm` | `double` | 残差ノルム `\|\|Ax - b\|\|` です。 |
| `R` | `MatrixObject?` | 列ピボット付き QR 分解で得られた上三角因子です。列はピボット後の順序で保存されます。`store_intermediates` が `true` の場合のみ保存され、それ以外では `null` です。 |
| `Qtb` | `double[]?` | 変換後の右辺ベクトル `Qᵀb` です。`store_intermediates` が `true` の場合のみ保存され、それ以外では `null` です。 |
| `Pivot` | `int[]?` | 列ピボット情報です。`Pivot[j]` は、ピボット後の第 `j` 列が元の何列目であったかを表します。`store_intermediates` が `true` の場合のみ保存され、それ以外では `null` です。 |
| `Tag` | `object?` | 追加のメタ情報を保持するための任意のタグ領域です。 |

#### 公開メソッド

| メソッド | 説明 |
| --- | --- |
| `ToString(bool omit = false)` | 結果オブジェクトの内容を読みやすい文字列に変換します。`omit` を `true` にすると、長い配列や行列の表示が省略されます。 |

#### 使用例

```csharp
var result = Solve(A, b, store_intermediates: true);

if (result.Status != LSQSolverStatus.Success)
{
    Console.WriteLine(result.Status);
    return;
}

double[] x = result.Solution;

Console.WriteLine($"Rank: {result.Rank}");
Console.WriteLine($"Residual norm: {result.ResidualNorm}");
Console.WriteLine(result.ToString(omit: true));
```

## ドキュメント

以下のドキュメントは今後追加予定です。

### 理論解説

[`theory_jp.md`](https://github.com/TaigaNakano/LSQSolver/blob/main/docs/theory_jp.md)

内容:

- 最小二乗問題
- 列ピボット付きQR分解
- 数値ランク判定
- 最小ノルム解
- Cholesky分解による解構成

### 応用例

[`polynomial-fit_jp.md`](https://github.com/TaigaNakano/LSQSolver/blob/main/docs/polynomial-fit_jp.md)

- 多項式近似
- Vandermonde行列
- 実用的なフィッティング例

[`gravity-inversion_jp.md`](https://github.com/TaigaNakano/LSQSolver/blob/main/docs/gravity-inversion_jp.md)

- 重力異常逆解析
- 劣決定最小二乗問題
- 最小ノルム解

### 性能評価

予備的なベンチマーク結果は以下にまとめています。

[`performance_jp.md`](https://github.com/TaigaNakano/LSQSolver/blob/main/docs/performance_jp.md)

このベンチマークでは、フルランクおよびランク落ちの密な正方最小二乗問題について、LSQSolver と GNU Octave を比較しています。

掲載している時間は10回実行したときの中央値であり、以下の環境で測定しています。

| 項目 | 値 |
|---|---|
| 機種名 | MacBook Pro |
| チップ | Apple M1 Pro |
| CPUコア数 | 合計8コア: パフォーマンスコア6、効率性コア2 |
| メモリ | 16 GB |

今回の測定では、以下の結果が得られました。

- フルランクの密行列では、このベンチマークにおいて `n >= 50` で LSQSolver は Octave の QR 分解より高速でした。
- `n = 2000` のフルランクケースでは、LSQSolver が `564.4 ms`、Octave の QR 分解が `1483.0 ms` でした。
- ランク落ち行列では、大きいサイズにおいて LSQSolver は Octave の QR 分解より高速でした。`n = 2000` では、LSQSolver が `764.7 ms`、Octave の QR 分解が `2102.2 ms` でした。
- Octave の `pinv` と比較すると、大きな行列では LSQSolver が大幅に高速でした。`n = 2000` では、ランク落ちケースで約 `12.8倍`、フルランクケースで約 `27.6倍` 高速でした。

これらの結果は予備的なものであり、ハードウェア、実行環境、コンパイラ設定、BLAS/LAPACK の構成、およびベンチマーク実装に依存します。

## 今後の対応予定

今後の開発は、ユーザーのニーズに応じて検討します。

候補としては、正則化付き最小二乗問題、複数右辺ベクトル、重み付き最小二乗問題への対応などがあります。これらは、ユーザー側での実装ミスを減らせる場合、不要なメモリ確保を避けられる場合、またはソルバー内部でアルゴリズム上の利点がある場合に対応を検討します。

## ライセンス

MIT License

