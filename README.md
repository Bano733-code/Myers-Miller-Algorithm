# Myers–Miller Linear-Space Sequence Alignment

A small C++ implementation of the **Myers–Miller algorithm for optimal sequence alignment in linear space**, based on the paper:

> Myers, E. W., & Miller, W. (1988). *Optimal alignments in linear space*. Bioinformatics, 4(1), 11–17.

## 📌 Project Overview

Sequence alignment is a fundamental problem in bioinformatics. Dynamic programming can find an optimal alignment, but the conventional approach requires a large `O(MN)` memory matrix for sequences of lengths `M` and `N`.

The Myers–Miller approach reduces the **space requirement to linear space** by combining:

* Gotoh's affine-gap alignment formulation
* Forward dynamic programming
* Reverse dynamic programming
* Divide-and-conquer
* Midpoint selection

The goal of this project is to understand and implement the main ideas behind the Myers–Miller algorithm in **C++**.

---

## 🧬 What the Algorithm Does

Given two sequences:

```text
Sequence A
Sequence B
```

the program finds an optimal alignment while using linear-space vectors for the main forward and reverse calculations.

The algorithm divides the alignment problem around a midpoint of sequence `A`.

The major steps are:

```text
Input sequences
       ↓
Forward calculation
       ↓
CC and DD vectors
       ↓
Reverse calculation
       ↓
RR and SS vectors
       ↓
Midpoint calculation
       ↓
Type 1 / Type 2 midpoint
       ↓
Divide into subproblems
       ↓
Recursive alignment
       ↓
Final optimal alignment
```

---

## 🧮 Scoring Model

This implementation uses an affine gap-cost model:

```text
gap(k) = g + h × k
```

where:

* `g` = gap opening cost
* `h` = gap extension cost
* `k` = gap length

The default parameters are:

```text
Match cost    = 0
Mismatch cost = 1
Gap opening   = 2
Gap extension = 0.5
```

These values follow the scoring setup used for demonstrating the algorithm in the paper.

---

## 🔬 Main Algorithm Components

### 1. Forward Calculation

The forward dynamic programming calculation produces:

```text
CC
DD
```

`CC[j]` represents the minimum alignment cost for the current prefix.

`DD[j]` represents the corresponding state ending in a deletion.

Only vectors are stored rather than the complete dynamic-programming matrix.

---

### 2. Reverse Calculation

The sequences are processed from the opposite direction to obtain:

```text
RR
SS
```

These vectors provide the cost of completing the alignment from the midpoint toward the end.

---

### 3. Midpoint Selection

For every possible position `j`, two possible midpoint types are evaluated.

#### Type 1

```text
CC[j] + RR[N-j]
```

#### Type 2

```text
DD[j] + SS[N-j] - g
```

The minimum value determines the midpoint used to divide the alignment problem.

The program displays these values so that the intermediate calculations can be inspected.

---

## 📊 Trace / Debug Output

The program can display the internal calculations:

```text
========== FORWARD ==========
CC: ...
DD: ...

========== REVERSE ==========
RR: ...
SS: ...

========== MIDPOINT ==========
j       Type 1       Type 2
0       ...          ...
1       ...          ...
2       ...          ...

Chosen midpoint:
i* = ...
j* = ...
Type = ...

========== FINAL ==========
Optimal alignment:
...
...

Minimum cost = ...
```

This trace is included specifically to make the algorithm easier to understand and verify.

---

## 📁 Project Structure

```text
MyersMiller/
│
├── main.cpp
│
├── scoring.h
├── scoring.cpp
│
├── forward.h
├── forward.cpp
│
├── reverse.h
├── reverse.cpp
│
├── midpoint.h
├── midpoint.cpp
│
├── alignment.h
└── alignment.cpp
```

### File Description

| File            | Purpose                                    |
| --------------- | ------------------------------------------ |
| `main.cpp`      | Takes input and runs the algorithm         |
| `scoring.cpp`   | Match, mismatch, and gap-cost functions    |
| `forward.cpp`   | Computes `CC` and `DD`                     |
| `reverse.cpp`   | Computes `RR` and `SS`                     |
| `midpoint.cpp`  | Evaluates Type 1 and Type 2 midpoint costs |
| `alignment.cpp` | Recursive divide-and-conquer alignment     |
| `.h` files      | Function declarations and data structures  |

---

## ⚙️ Requirements

* C++17 or later
* GNU `g++` compiler

No external libraries are required.

The project can also be compiled in **Google Colab** using its built-in C++ compiler.

---

## ▶️ Compilation

Compile all source files together:

```bash
g++ -std=c++17 -O2 main.cpp scoring.cpp forward.cpp reverse.cpp midpoint.cpp alignment.cpp -o myers_miller
```

Run:

```bash
./myers_miller
```

---

## 💻 Example Input

```text
Enter sequence A: ACGT
Enter sequence B: AGT
```

The program then displays the forward vectors, reverse vectors, midpoint calculations, and final alignment.

---

## 🎯 Educational Purpose

This project was developed as a **small algorithm implementation project** to understand the Myers–Miller linear-space sequence alignment method.

The focus is not on building a large application, but on understanding how:

1. affine-gap dynamic programming works,
2. linear-space vectors replace the full DP matrix,
3. forward and reverse calculations identify a midpoint,
4. the problem is divided recursively, and
5. an optimal alignment is reconstructed.

---

## 📚 Reference

Myers, E. W., & Miller, W. (1988).

**Optimal alignments in linear space.**

*Bioinformatics*, 4(1), 11–17.

---

## ⚠️ Implementation Note

This repository is an **educational implementation inspired by the Myers–Miller method**.

The implementation demonstrates the core linear-space forward/reverse calculations, midpoint selection, and divide-and-conquer strategy. Small terminal subproblems are solved using a full dynamic-programming procedure for simplicity.

It should therefore not be interpreted as a line-by-line reproduction of the original software implementation described in the paper.
