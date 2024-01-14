pub mod circulant;
pub mod toeplitz;
pub mod util;

use circulant::*;
use util::*;

pub trait Matrix<F> {
    fn get_len(&self) -> usize;
    fn get_matrix(&self) -> Vec<F>;
}
