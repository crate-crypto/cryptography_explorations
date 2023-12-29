use ark_bn254::{Fr, G1Projective as G1};
use ark_ec::Group;
use ark_ff::PrimeField;
use ark_std::UniformRand;
use ark_std::Zero;
use rand::rngs::ThreadRng;
use std::error::Error;
use std::ops::Add;
use std::ops::Mul;

// NOTE References:
///**
//  * Create a pedersen commitment (point) from an array of input fields.
//  * Left pads any inputs less than 32 bytes.
//  */
//
// The input is Vec<Vec<u8>> in Rust.
// export function pedersenCommit(input: Buffer[]) {
//     if (!input.every(i => i.length <= 32)) {
//       throw new Error('All input buffers must be <= 32 bytes.');
//     }
//
// Left padding:
//     input = input.map(i => (i.length < 32 ? Buffer.concat([Buffer.alloc(32 - i.length, 0), i]) : i));
//     const point = BarretenbergSync.getSingleton().pedersenCommit(input.map(i => new Fr(i)));
// toBuffer returns Uint8Arrays (browser/worker-boundary friendly).
//     return [Buffer.from(point.x.toBuffer()), Buffer.from(point.y.toBuffer())];
//   }
//
// Implementation itself:
// typename Curve::AffineElement pedersen_commitment_base<Curve>::commit_native(const std::vector<Fq>& inputs,
// const GeneratorContext context)
// {
//     const auto generators = context.generators->get(inputs.size(), context.offset, context.domain_separator);
//     Element result = Group::point_at_infinity;
//
//     for (size_t i = 0; i < inputs.size(); ++i) {
//         result += Element(generators[i]) * static_cast<uint256_t>(inputs[i]);
//     }
//     return result.normalize();
// }

// Derive generators:
// inline static constexpr std::vector<affine_element> derive_generators(
//     const std::vector<uint8_t>& domain_separator_bytes,
//     const size_t num_generators,
//     const size_t starting_index = 0)
// {
//     std::vector<affine_element> result;
//     const auto domain_hash = blake3::blake3s_constexpr(&domain_separator_bytes[0], domain_separator_bytes.size());
//     std::vector<uint8_t> generator_preimage;
//     generator_preimage.reserve(64);
//     std::copy(domain_hash.begin(), domain_hash.end(), std::back_inserter(generator_preimage));
//     for (size_t i = 0; i < 32; ++i) {
//         generator_preimage.emplace_back(0);
//     }
//     for (size_t i = starting_index; i < starting_index + num_generators; ++i) {
//         auto generator_index = static_cast<uint32_t>(i);
//         uint32_t mask = 0xff;
//         generator_preimage[32] = static_cast<uint8_t>(generator_index >> 24);
//         generator_preimage[33] = static_cast<uint8_t>((generator_index >> 16) & mask);
//         generator_preimage[34] = static_cast<uint8_t>((generator_index >> 8) & mask);
//         generator_preimage[35] = static_cast<uint8_t>(generator_index & mask);
//         result.push_back(affine_element::hash_to_curve(generator_preimage));
//     }
//     return result;
// }

const N: usize = 32/* replace with your desired size */;

#[derive(Clone)]
pub struct GeneratorContext {
    pub generators: Vec<Vec<G1>>, // Replace it with generics
}

use std::marker::PhantomData;
pub struct Commitment<C: Group> {
    group: PhantomData<C>,
}

// TODO: Think on how to make it faster if possible
fn left_pad_bytes(input: &[u8], length: usize, padding_byte: u8) -> Vec<u8> {
    let padding_len = length.saturating_sub(input.len());
    let mut padded = vec![padding_byte; padding_len];
    padded.extend_from_slice(input);
    padded
}

fn transform_to_array_slice(input: &[Vec<u8>]) -> Result<Vec<[u8; N]>, Box<dyn Error>> {
    input
        .iter()
        .map(|vec| {
            if vec.len() > 32 {
                Err("Inner vector length should be less than or equal to 32")?
            } else {
                let padded_bytes = left_pad_bytes(vec.as_slice(), N, 0);
                let mut array = [0u8; N];
                array.copy_from_slice(&padded_bytes[..N]);
                Ok(array)
            }
        })
        .collect()
}

fn transform_to_fr(input: Vec<[u8; N]>) -> Vec<Fr> {
    let mut result = Vec::new();
    for item in input {
        result.push(Fr::from_le_bytes_mod_order(&item));
    }
    result
}

// TODO: use pregenerated generators rather than classic C = sG + rH
fn pedersen_commit_based_on_fr(input: &[Fr], blinding_factor: Fr) -> G1 {
    let commitment_point = input
        .iter()
        .map(|fi| G1::generator().mul(fi))
        .fold(G1::zero(), |acc, point| acc.add(&point))
        .add(&G1::generator().mul(&blinding_factor));

    commitment_point
}

pub fn pedersen_commit(input: &[Vec<u8>], mut rng: ThreadRng) -> Result<G1, Box<dyn Error>> {
    // Random value for blinding factor
    let r: Fr = Fr::rand(&mut rng);

    if input.len().is_zero() {
        return Err("The array is empty".into());
    }

    let input_padded: Vec<[u8; N]> = transform_to_array_slice(input).unwrap();

    let input_fr = transform_to_fr(input_padded);

    // Pedersen commitment
    let commitment_point = pedersen_commit_based_on_fr(&input_fr, r);
    Ok(commitment_point)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_left_pad_bytes() {
        // Padding needed
        let input1 = [1, 2, 3];
        let padded_result1 = left_pad_bytes(&input1, 6, 0);
        assert_eq!(padded_result1, vec![0, 0, 0, 1, 2, 3]);

        // No padding needed
        let input2 = [4, 5, 6];
        let padded_result2 = left_pad_bytes(&input2, 3, 0);
        assert_eq!(padded_result2, vec![4, 5, 6]);

        // Empty input
        let input3 = [];
        let padded_result3 = left_pad_bytes(&input3, 4, 9);
        assert_eq!(padded_result3, vec![9, 9, 9, 9]);

        // Zero length
        let input4 = [7, 8, 9];
        let padded_result4 = left_pad_bytes(&input4, 0, 3);
        assert_eq!(padded_result4, vec![7, 8, 9]);

        // Length equals input length
        let input5 = [1, 2, 3];
        let padded_result5 = left_pad_bytes(&input5, 3, 9);
        assert_eq!(padded_result5, vec![1, 2, 3]);
    }

    #[test]
    // Non-empty commit
    fn test_pedersen_commit() {
        // Non-empty input
        let input1 = vec![vec![1, 2, 3], vec![4, 5, 6], vec![7, 8, 9]];
        let rng1 = rand::thread_rng();
        let commitment1 = pedersen_commit(&input1, rng1);
        assert!(!commitment1.unwrap().is_zero());

        // Empty input
        let input2 = vec![];
        let rng2 = rand::thread_rng();
        let commitment2 = pedersen_commit(&input2, rng2);
        assert!(commitment2.is_err());

        // Random input
        let input3 = vec![vec![10, 20, 30], vec![40, 50, 60]];
        let rng3 = rand::thread_rng();
        let commitment3 = pedersen_commit(&input3, rng3);
        assert!(!commitment3.unwrap().is_zero());
    }
}
