# Number Theoretic Transform (NTT) in Rust

Dieses Projekt implementiert die **Number Theoretic Transform (NTT)** und deren inverse Transformation in Rust. Die NTT ist eine Variante der diskreten Fourier-Transformation (DFT), die auf modularer Arithmetik basiert und häufig in der Kryptographie sowie schnellen Polynom-Multiplikation eingesetzt wird.

---

## Funktionen

- **Forward NTT (rekursiv)**  
  Berechnet die Vorwärts-NTT eines Eingabevektors mit einem primitiven n-ten Einheitswurzel modulo eines Primzahlmodulus.

- **Inverse NTT (rekursiv)**  
  Berechnet die inverse Transformation mittels modularer Inversion.

- **Primitive n-te Einheitswurzel finden**  
  Findet eine primitive n-te Einheitswurzel modulo des gegebenen Modulus.

- **Hilfsfunktionen**  
  Modulare Potenzierung, modulare Multiplikation, erweiterter euklidischer Algorithmus, etc.

---

## Voraussetzungen

- Rust (empfohlen: aktuelle stabile Version)
- Git (zum Klonen des Repositories)

---

## Installation und Nutzung

1. Repository klonen:

   ```bash
   git clone https://github.com/HasanAydin004/NumberTheoreticTransform.git
   cd NumberTheoreticTransform
2. Projekt bauen:

   ```bash
   cargo build --release
3. Benchmarks ausführen:
   ```bash
   cargo bench
Beispiel Code : 
```rust
let input = vec![1, 2, 3, 4, 0, 0, 0, 0];
let modulus = 97;
let totient = 96;
let root = find_primitive_nth_root(input.len() as i128, totient, modulus);

let ntt_result = forward_ntt_recursive(&input, root, modulus);
let inverse_result = inverse_ntt_recursive(&ntt_result, root, modulus);
