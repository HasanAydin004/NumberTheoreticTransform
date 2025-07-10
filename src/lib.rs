// Computes the forward Number-Theoretic Transform (NTT) of the input vector `invec`,
// using the provided primitive n-th root of unity `root`, modulo `modulus`.
pub fn forward_ntt_recursive(invec: &Vec<i128>, root: i128, modulus: i128) -> Vec<i128> {
    let n = invec.len();
    if n == 1 {
        return vec![invec[0] % modulus];
    }
    let half = n / 2;

    // Split into even and odd indices
    let even: Vec<i128> = invec.iter().step_by(2).cloned().collect();
    let odd: Vec<i128> = invec.iter().skip(1).step_by(2).cloned().collect();

    // Square the root to pass down recursively
    let root_squared = modular_pow(root, 2, modulus);
    let even_ntt = forward_ntt_recursive(&even, root_squared, modulus);
    let odd_ntt = forward_ntt_recursive(&odd, root_squared, modulus);

    // Combine step
    let mut out = vec![0; n];
    let mut omega = 1;

    for i in 0..half {
        let t = (omega * odd_ntt[i]) % modulus;
        out[i] = (even_ntt[i] + t + modulus) % modulus;
        out[i + half] = (even_ntt[i] - t + modulus) % modulus;
        omega = (omega * root) % modulus;
    }

    out
}

// Computes the inverse Number-Theoretic Transform (INTT) of the input vector `invec`,
// using the provided primitive n-th root of unity `root`, modulo `modulus`.
//
// The inverse is computed by applying a forward NTT with the reciprocal root
// and then multiplying each coefficient by the modular inverse of the vector length.
pub fn inverse_ntt_recursive(invec: &Vec<i128>, root: i128, modulus: i128) -> Vec<i128> {
    let n = invec.len() as i128;
    let inv_root = modular_inverse(root, modulus);
    let mut outvec = forward_ntt_recursive(invec, inv_root, modulus);

    let scaler = modular_inverse(n, modulus);
    for i in 0..outvec.len() {
        outvec[i] = modular_mul(outvec[i], scaler, modulus);
    }

    outvec
}




// Returns a primitive `degree`-th root of unity modulo `modulus`.
// Requires that `totient` (Euler's totient of modulus) is divisible by `degree`.
// If `modulus` is prime, such a root is guaranteed to exist.
pub fn find_primitive_nth_root(degree: i128, totient: i128, modulus: i128) -> i128 {
    let generator = find_multiplicative_generator(totient, modulus);
    modular_pow(generator, totient / degree, modulus)
}

// Finds a generator (primitive root) of the multiplicative group modulo `modulus`.
// Assumes `totient` = φ(modulus). For prime `modulus`, such a generator always exists.
fn find_multiplicative_generator(totient: i128, modulus: i128) -> i128 {
    for candidate in 1..modulus {
        if is_primitive_nth_root(candidate, totient, modulus) {
            return candidate;
        }
    }
    panic!("No primitive root found. Are you sure the modulus is correct?");
}

// Checks whether `val` is a primitive root modulo `modulus`.
// That means: val^totient ≡ 1 (mod modulus), but not for any smaller exponent dividing totient.
fn is_primitive_nth_root(val: i128, totient: i128, modulus: i128) -> bool {
    // Must satisfy Fermat's little theorem
    if modular_pow(val, totient, modulus) != 1 {
        return false;
    }

    // Must NOT satisfy for any smaller divisor totient/p (p prime)
    for p in get_unique_prime_factors(totient) {
        if modular_pow(val, totient / p, modulus) == 1 {
            return false;
        }
    }

    true
}

// Returns the unique prime factors of `n` in ascending order.
// Example: unique_prime_factors(60) = [2, 3, 5]
fn get_unique_prime_factors(mut n: i128) -> Vec<i128> {
    let mut result = Vec::new();
    let mut i = 2;

    while i * i <= n {
        if n % i == 0 {
            result.push(i);
            while n % i == 0 {
                n /= i;
            }
        }
        i += 1;
    }

    if n > 1 {
        result.push(n);
    }

    result
}
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

fn modular_mul(lhs: i128, rhs: i128, modulo: i128) -> i128 {
    ((lhs % modulo) * (rhs % modulo)) % modulo
}

// Computes the modular inverse of `a` modulo `modulus`, i.e., x such that (a * x) % modulus == 1.
//Panics if no inverse exists.
pub fn modular_inverse(a: i128, modulus: i128) -> i128 {
    let (g, x, _) = compute_extended_gcd(a, modulus);
    if g != 1 {
        panic!("No modular inverse exists for {} mod {}", a, modulus);
    }
    (x % modulus + modulus) % modulus
}
//Extended Euclidean Algorithm: returns gcd(a, b) and Bézout coefficients (x, y)
fn compute_extended_gcd(a: i128, b: i128) -> (i128, i128, i128) {
    if b == 0 {
        (a, 1, 0)
    } else {
        let (g, x1, y1) = compute_extended_gcd(b, a % b);
        (g, y1, x1 - (a / b) * y1)
    }
}
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
pub fn two_poly(ntt_a: &[i128], ntt_b: &[i128], modulus: i128) -> Vec<i128> {
    assert_eq!(ntt_a.len(), ntt_b.len(), "Polynome müssen gleich lang sein");

    ntt_a.iter()
        .zip(ntt_b.iter())
        .map(|(&x, &y)| modular_mul(x, y, modulus))
        .collect()
}