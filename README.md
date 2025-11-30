# Visual Grouping Tool (Efficient Hierarchical Clustering)

**Assignment 03** **Author:** Muhammad Obaidullah [24L-0509]  
**Language:** C++  

---

## 📌 Project Overview

This project implements **Agglomerative Hierarchical Clustering (Single Linkage)** for 2D spatial data. 

Standard implementations of this algorithm often run in **$O(N^3)$** time because they repeatedly scan a distance matrix to find the closest pair. This project optimizes the process to **$O(N^2 \log N)$** by utilizing two key data structures:
1.  **Min-Heap:** For $O(1)$ access to the closest pair of groups.
2.  **Lookup Matrix:** For $O(1)$ access to the location of specific edges inside the Heap during updates.

## 🚀 Key Features

* **Custom Min-Heap:** Implemented from scratch (without STL `priority_queue`) to support `updateKey` operations.
* **Position Lookup Matrix:** A 2D mapping (`matrix[u][v] = index`) that tracks exactly where the distance between any two groups is stored in the heap array.
* **Efficient Merging:** Uses `std::list::splice` to merge group members in $O(1)$ constant time without copying data.
* **Lazy Deletion Strategy:** Handles merged groups by marking them inactive rather than resizing the matrix.

## ⚙️ How It Works (The Logic)

1.  **Initialization:** * Reads $N$ points.
    * Computes $N(N-1)/2$ pairwise distances.
    * Builds a **Min-Heap** of these distances.
    * Populates the **Lookup Matrix**.

2.  **Clustering Loop:** * **Extract:** The heap returns the closest pair of groups $(U, V)$.
    * **Merge:** Points from Group $V$ are spliced into Group $U$. Group $V$ is deactivated.
    * **Update:** For every other active group $K$:
        * Calculate new distance: $Min(Dist(U,K), Dist(V,K))$.
        * Use **Lookup Matrix** to find the specific node for $(U, K)$ in the heap.
        * **Bubble Up/Down:** Update the priority in the heap.
        * "Delete" the old $(V, K)$ connection by setting it to Infinity.

## 📄 Input Format

The program expects a file named `inputData.txt` in the same directory.

**Format:**
```text
<N_Points> , <M_Target_Groups>
<X1> <Y1>
<X2> <Y2>
...
