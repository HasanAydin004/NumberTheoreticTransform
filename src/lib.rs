/// Führt die Vorwärts-Number-Theoretic Transform (NTT) des Vektors `invec` durch.
///
/// Die Implementierung basiert auf dem rekursiven Radix-2–Cooley–Tukey-Algorithmus
/// (siehe Abschnitt 3.3 im Paper). Diese Methode reduziert die Komplexität
/// der Polynomtransformation von O(n²) auf O(n log n).
///
/// # Parameter:
/// - `invec`: Eingangsvektor (Polynomkoeffizienten a₀, a₁, ..., aₙ₋₁)
/// - `root`: Primitive n-te Einheitswurzel ω ∈ Z_q (siehe Abschnitt 2.2)
/// - `modulus`: Primzahl q > n, für Modulo-Arithmetik in Z_q
///
/// # Rückgabe:
/// Vektor der NTT-Koeffizienten Â₀, Â₁, ..., Âₙ₋₁ im Frequenzbereich.
///
/// # Mathematische Grundlage:
/// Die NTT entspricht der diskreten Fourier-Transformation über dem endlichen Ring Z_q:
/// Âⱼ = Σᵢ aᵢ · ω^(i·j) mod q  (vgl. Abschnitt 3.2)
pub fn forward_ntt_recursive(invec: &Vec<i128>, root: i128, modulus: i128) -> Vec<i128> {
    let n = invec.len();

    // Rekursionsbasis: Ein einzelnes Element ⇒ einfach modulo q nehmen
    if n == 1 {
        return vec![invec[0] % modulus];
    }

    let half = n / 2;

    // Zerlege den Vektor in gerad- und ungerad-indizierte Elemente
    // Entspricht der Teilung in zwei kleinere NTTs (Divide-and-Conquer)
    let even: Vec<i128> = invec.iter().step_by(2).cloned().collect();  // a₀, a₂, a₄, ...
    let odd: Vec<i128> = invec.iter().skip(1).step_by(2).cloned().collect();  // a₁, a₃, a₅, ...

    // Für die rekursive Unter-NTT muss ω durch ω² ersetzt werden
    // (siehe Cooley–Tukey-Struktur: ω² hat Ordnung n/2)
    let root_squared = modular_pow(root, 2, modulus);

    // Rekursiver Aufruf für die Teilstücke
    let even_ntt = forward_ntt_recursive(&even, root_squared, modulus);
    let odd_ntt = forward_ntt_recursive(&odd, root_squared, modulus);

    // Initialisiere Ausgabevektor mit Nullen
    let mut out = vec![0; n];
    let mut omega = 1;  // ω⁰ = 1, Startwert

    // Kombinationsphase ("Butterfly"-Operation, siehe Abschnitt 3.3):
    // Berechne die NTT-Werte aus den beiden halben NTTs
    for i in 0..half {
        // Berechne t = ω^i · odd_ntt[i] mod q
        // (ω wird schrittweise auf ω¹, ω², ..., ω^(n/2 - 1) erhöht)
        let t = (omega * odd_ntt[i]) % modulus;

        // Obere Butterfly-Zweig: even + t mod q
        out[i] = (even_ntt[i] + t + modulus) % modulus;

        // Untere Butterfly-Zweig: even - t mod q
        // (+modulus vermeidet negative Zwischenwerte)
        out[i + half] = (even_ntt[i] - t + modulus) % modulus;

        // Aktualisiere ω: ω ← ω · root mod q
        // (das entspricht dem nächsten Potenzwert ωⁱ⁺¹)
        omega = (omega * root) % modulus;
    }

    // Transformierter Vektor im Frequenzbereich (NTT)
    out
}



/// Führt die **inverse Number-Theoretic Transform (INTT)** auf dem Eingabevektor `invec` aus.
///
/// Die Inverse-NTT berechnet das ursprüngliche Polynom zurück aus seinem Frequenzbereich.
/// Dies geschieht in zwei Schritten:
///
/// 1. Durchführung einer Vorwärts-NTT mit der **inversen** primitiven Einheitswurzel `ω⁻¹`,
/// 2. anschließende **Skalierung** jedes Koeffizienten mit dem multiplikativen Inversen von `n` modulo `q`,
///    um die korrekte Rücktransformation sicherzustellen.
///
/// Mathematischer Hintergrund (siehe Abschnitt 3.2 des Papers):
/// Die INTT ist definiert als:
///     aᵢ = (1/n) * Σⱼ=0ⁿ⁻¹ Âⱼ · ω⁻ⁱʲ  mod q
///
/// - `invec`: Transformierter Vektor Â im Frequenzbereich (z. B. Ergebnis der Vorwärts-NTT)
/// - `root`: Die ursprüngliche primitive n-te Einheitswurzel ω ∈ Z_q
/// - `modulus`: Das verwendete q (typischerweise eine große Primzahl)
///
/// Gibt den Zeitbereichsvektor a zurück.
pub fn inverse_ntt_recursive(invec: &Vec<i128>, root: i128, modulus: i128) -> Vec<i128> {
    let n = invec.len() as i128;

    // Berechne das multiplikative Inverse der Einheitswurzel ω ⇒ ω⁻¹ ∈ Z_q
    let inv_root = modular_inverse(root, modulus);

    // Berechne NTT mit inverser Wurzel ⇒ entspricht der unskalierten INTT
    let mut outvec = forward_ntt_recursive(invec, inv_root, modulus);

    // Normierungsschritt: Skaliere jedes Ergebnis mit dem Inversen von n modulo q ⇒ 1/n mod q
    let scaler = modular_inverse(n, modulus);

    // Jeder Koeffizient wird mit (1/n mod q) multipliziert, um korrekte Rücktransformation zu erhalten
    for i in 0..outvec.len() {
        outvec[i] = modular_mul(outvec[i], scaler, modulus);
    }

    outvec
}





/// Findet eine primitive n-te Einheitswurzel ω ∈ Z_q.
/// Voraussetzung: n | φ(q), also n teilt q - 1.
///
/// Entspricht Abschnitt 2.2 & 4.1 im Paper.
pub fn find_primitive_nth_root(degree: i128, totient: i128, modulus: i128) -> i128 {
    let generator = find_multiplicative_generator(totient, modulus);
    modular_pow(generator, totient / degree, modulus)
}


/// Sucht einen Generator g ∈ Z_q*, also ein Element der Ordnung φ(q) = q-1.
/// Verwendet in der Einheitswurzel-Suche (Abschnitt 2.3).
fn find_multiplicative_generator(totient: i128, modulus: i128) -> i128 {
    for candidate in 1..modulus {
        if is_primitive_nth_root(candidate, totient, modulus) {
            return candidate;
        }
    }
    panic!("Kein primitiver Generator gefunden – ist q wirklich eine Primzahl?");
}

/// Prüft, ob `val` ein primitiver Generator modulo `modulus` ist,
/// d. h. ein Erzeuger der multiplikativen Gruppe Z*_q.
///
/// Hintergrund:
/// Ein Element `g` ∈ Z*_q ist genau dann ein primitiver Generator, wenn es Ordnung φ(q) = q−1 hat.
/// Das bedeutet: Seine Potenzen durchlaufen alle invertierbaren Reste modulo q.
///
/// Dafür müssen zwei Bedingungen erfüllt sein:
/// 1. g^φ(q) ≡ 1 mod q       (Fermats kleiner Satz)
/// 2. Für alle Primfaktoren p von φ(q): g^(φ(q)/p) ≠ 1 mod q
///
/// Diese Funktion wird bei der Suche nach einer primitiven Einheitswurzel verwendet (vgl. Abschnitt 2.3).
fn is_primitive_nth_root(val: i128, totient: i128, modulus: i128) -> bool {
    // 1. Prüfe g^φ(q) ≡ 1 mod q
    // Dies muss immer gelten – folgt direkt aus Fermats kleinem Satz.
    if modular_pow(val, totient, modulus) != 1 {
        return false;
    }

    // 2. Prüfe, ob val keine kleinere Ordnung hat:
    // Für alle Primfaktoren p von φ(q) gilt:
    // Wenn val^(φ(q)/p) ≡ 1 (mod q), dann liegt val in einer echten Untergruppe ⇒ kein Generator.
    for p in get_unique_prime_factors(totient) {
        if modular_pow(val, totient / p, modulus) == 1 {
            return false; // val liegt in einer echten Untergruppe ⇒ disqualifiziert
        }
    }

    // Wenn alle Tests bestanden sind, ist val ein Generator von Z*_q
    true
}



/// Gibt die eindeutigen Primfaktoren von `n` in aufsteigender Reihenfolge zurück.
/// Wird z.B. in `is_primitive_nth_root` verwendet, um zu prüfen, ob ein Kandidat
/// ein echter Generator (primitive Wurzel) ist.
///
/// Beispiel: get_unique_prime_factors(60) = [2, 3, 5]
///
/// Mathematischer Kontext:
/// Um zu prüfen, ob ein Element g ∈ Z*_q ein Generator ist, muss man alle Primfaktoren
/// von φ(q) = q − 1 kennen, um g^(φ(q)/p) ≠ 1 mod q zu prüfen (siehe Abschnitt 2.3).
fn get_unique_prime_factors(mut n: i128) -> Vec<i128> {
    let mut result = Vec::new(); // Enthält die gefundenen Primfaktoren ohne Wiederholungen
    let mut i = 2; // Beginne mit der kleinsten Primzahl

    // Solange i * i ≤ n, d.h. wir testen nur bis zur Quadratwurzel von n
    while i * i <= n {
        // Wenn i ein Teiler von n ist, ist i ein Primfaktor
        if n % i == 0 {
            result.push(i); // i zur Ergebnisliste hinzufügen

            // Entferne alle Vielfachen von i aus n (z. B. 60 → 30 → 15 → 5 für i = 2)
            while n % i == 0 {
                n /= i;
            }
        }
        i += 1; // Nächster Kandidat
    }

    // Falls nach der Schleife n > 1 ist, bleibt ein letzter Primfaktor übrig (z. B. eine große Primzahl)
    if n > 1 {
        result.push(n);
    }

    result
}

// Schnelle Exponentiation: base^exp mod modulus (log Zeit)
fn modular_pow(base: i128, exp: i128, modulus: i128) -> i128 {
    let mut result = 1;
    let mut base = base % modulus;//mann darf vorher schon die basis kürzen damit man später leichter rechnen kann
    let mut exp = exp;
    while exp > 0 {
        if exp % 2 == 1 {
            result = result * base % modulus;
        }
        base = base * base % modulus;
        exp /= 2;
    }
    result
}

// Multiplikation zweier Zahlen modulo q
fn modular_mul(lhs: i128, rhs: i128, modulo: i128) -> i128 {
    ((lhs % modulo) * (rhs % modulo)) % modulo
}
/// Modularinverse: Findet x mit (a * x) ≡ 1 mod modulus
pub fn modular_inverse(a: i128, modulus: i128) -> i128 {
    let (g, x, _) = compute_extended_gcd(a, modulus);
    if g != 1 {
        panic!("Kein modulares Inverses für {} mod {}", a, modulus);
    }
    (x % modulus + modulus) % modulus
}
// Erweiterter Euklid: gcd(a, b) und Bézout-Koeffizienten
fn compute_extended_gcd(a: i128, b: i128) -> (i128, i128, i128) {
    if b == 0 {
        (a, 1, 0)
    } else {
        let (g, x1, y1) = compute_extended_gcd(b, a % b);
        (g, y1, x1 - (a / b) * y1)
    }
}
/// Klassische, langsame NTT durch direkte Matrixmultiplikation:
/// Siehe Gleichung in Abschnitt 3.2 – O(n²)-Algorithmus.
pub fn ntt_based_on_omega(a: &[i128], omega: i128, q: i128) -> Vec<i128> {
    let n = a.len();
    let mut result = vec![0i128; n];
    for j in 0..n {
        let mut sum = 0i128;
        for i in 0..n {
            let power = modular_pow(omega, (i * j) as i128, q);
            sum = (sum + a[i] * power) % q;
        }
        result[j] = sum;
    }
    result
}
/// Multipliziert zwei Polynome im Frequenzbereich (komponentenweise).
/// Verwendet nach Vorwärts-NTT: Ĉ = Â ◦ B̂.
/// Siehe Abschnitt 4.5.
pub fn two_poly(ntt_a: &[i128], ntt_b: &[i128], modulus: i128) -> Vec<i128> {
    assert_eq!(ntt_a.len(), ntt_b.len(), "Polynome müssen gleich lang sein");

    ntt_a.iter()
        .zip(ntt_b.iter())
        .map(|(&x, &y)| modular_mul(x, y, modulus))
        .collect()
}
