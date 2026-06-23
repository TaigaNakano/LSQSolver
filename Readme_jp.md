# LSQSolver

.NET向け密行列最小二乗ソルバーです。

以下の問題をサポートします。

- 過決定問題
- 劣決定問題
- ランク落ち問題
- 最小ノルム解

SVDを使用せず、列ピボット付きQR分解（CPQR）に基づいて解を求めます。

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

### 戻り値

`Solve()` は `LSQSolverResult` を返します。

`result.Solution` を使う前に `result.Status` を確認してください。

```csharp
if (result.Status == LSQSolverStatus.Success)
{
    double[] x = result.Solution;
}
```

### ステータス確認

```csharp
if (result.Status != LSQSolverStatus.Success)
{
    Console.WriteLine(result.Status);
    return;
}
```

#### ステータスコード

想定されるステータス値:

- `Success`
- `NullMatrix`
- `EmptyMatrix`
- `NullVector`
- `DimensionMismatch`
- `InvalidMatrix`
- `InvalidVector`
- `CholeskyFailed`


### 解と診断情報の取得

```csharp
double[] x = result.Solution;

Console.WriteLine(result.Rank);
Console.WriteLine(result.ResidualNorm);
```

## ドキュメント

以下のドキュメントは今後追加予定です。

### 理論解説

`docs/theory.md`

内容:

- 最小二乗問題
- 列ピボット付きQR分解
- 数値ランク判定
- 最小ノルム解
- Cholesky分解による解構成

### 応用例

`docs/polynomial-fit.md`

- 多項式近似
- Vandermonde行列
- 実用的なフィッティング例

`docs/gravity-inversion.md`

- 重力異常逆解析
- 劣決定最小二乗問題
- 最小ノルム解

### 性能評価

性能評価

予備的なベンチマーク結果は以下にまとめています。

docs/performance_jp.md

このベンチマークでは、フルランクおよびランク落ちの密な正方最小二乗問題について、LSQSolver と GNU Octave を比較しています。

掲載している時間は10回実行したときの中央値であり、以下の環境で測定しています。

| 項目     | 値                        |
| ------ | ------------------------ |
| 機種名    | MacBook Pro              |
| チップ    | Apple M1 Pro             |
| CPUコア数 | 合計8コア: パフォーマンスコア6、効率性コア2 |
| メモリ    | 16 GB                    |

今回の測定では、以下の結果が得られました。

フルランクの密行列では、n >= 50 において LSQSolver は Octave の QR 分解より高速でした。
n = 2000 のフルランクケースでは、LSQSolver が 1134.0 ms、Octave の QR 分解が 1483.0 ms でした。
ランク落ち行列では、小さいサイズでは Octave の QR 分解が高速でしたが、大きいサイズでは LSQSolver が高速になりました。n = 2000 では、LSQSolver が 1701.6 ms、Octave の QR 分解が 2102.2 ms でした。
Octave の pinv と比較すると、大きな行列では LSQSolver が大幅に高速でした。n = 2000 では、ランク落ちケースで約 5.7倍、フルランクケースで約 13.8倍 高速でした。

これらの結果は予備的なものであり、ハードウェア、実行環境、コンパイラ設定、BLAS/LAPACK の構成、およびベンチマーク実装に依存します。

状況: 作成中

## ライセンス

MIT License

