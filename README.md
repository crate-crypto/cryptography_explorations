# Cryptography exploration

The main goal of this project is to implement fast KZG commitments and test other crypto.

## Lagrange interpolation

This part of the project is a Lagrange interpolation algorithm adaptation for use with points on BN254 curve. The primary goal is to test opportunities of ark-poly crate, that doesn't include constructing polynomials based on points. Lagrange interpolation is a mathematical technique that finds a polynomial passing through given points. Formula:

![Lagrange Interpolation](https://pbs.twimg.com/media/D13nksVW0AMo9EU.jpg)
