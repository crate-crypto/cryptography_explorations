// let p_point = s_poly.evaluate(&Fr::from(3));
// println!("p_point {:?}", p_point);

//
// let p_poly_un = DenseOrSparsePolynomial::from(p_poly);
// let q_poly_un = DenseOrSparsePolynomial::from(q_poly);
// let (quot, rem) = p_poly_un.divide_with_q_and_r(&q_poly_un).unwrap();
// println!("{:?}", quot);

// assert_eq!(s_poly, quot);
// assert!(rem.is_empty());

// // Compute pairings
// let pairing1 = Bls12_381::pairing(p_point * g1_generator, q_inv_point);

// // Perform polynomial division
// let division_poly = &p_poly / &q_poly;
// let div_point = evaluate_polynomial_on_curve(&division_poly, &g1_generator);

// // Compute pairings for the divided polynomial
// let pairing2 = Bls12_381::pairing(div_point, g2_generator);

// // Assert that the pairings are equal
// assert_eq!(pairing1, pairing2);

// // Random G1 point
// let p_scalar = Fr::rand(&mut rng);
// println!("{}", p_scalar);
// let p = g1_generator * &p_scalar;

// // Random G2 point
// let q_scalar = Fr::rand(&mut rng);
// println!("{}", q_scalar);

// // Invert the scalar
// let q_scalar_inv = q_scalar.inverse().unwrap();
// let q_inv = g2_generator * &q_scalar_inv;

// // Multiply
// let pairing1 = Bls12_381::pairing(p, q_inv);

// let division = p_scalar / q_scalar;
// println!("{}", division);
// let div_point = g1_generator * &division; // Возмущается, если писать не в том порядке
// let pairing2 = Bls12_381::pairing(div_point, g2_generator);
// assert_eq!(pairing1, pairing2);
