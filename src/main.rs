use number_theoretic_transform::{forward_ntt_recursive, inverse_ntt_recursive};

fn main() {
    let input: Vec<i128> = (0..128).collect();
    let root = 17;
    let modulus = 3329;
    let ntt = forward_ntt_recursive(&input, root, modulus);
    println!("{:?}", ntt);
    let invntt = inverse_ntt_recursive(&ntt, root, modulus);
    if input == invntt {
        println!("Alles passt");
        println!("{:?}", inverse_ntt_recursive(&ntt,root, modulus));
    }
}