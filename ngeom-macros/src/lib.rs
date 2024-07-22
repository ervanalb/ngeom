extern crate proc_macro;
use nalgebra::DVector;
use proc_macro2::Span;
use proc_macro2::TokenStream;
use quote::quote;
use syn::parse::{Parse, ParseStream, Result};
use syn::punctuated::Punctuated;
use syn::{parse_macro_input, Ident, LitInt, Token};

struct Input {
    basis_vector_squares: Vec<isize>,
}

impl Parse for Input {
    fn parse(input: ParseStream) -> Result<Self> {
        let lit_list = Punctuated::<LitInt, Token![,]>::parse_terminated(input)?;
        Ok(Input {
            basis_vector_squares: lit_list
                .iter()
                .map(|lit| lit.base10_parse::<isize>())
                .collect::<Result<Vec<_>>>()?,
        })
    }
}

fn gen_algebra2(input: Input) -> TokenStream {
    let Input {
        basis_vector_squares,
    } = input;

    fn bubble_sort_count_swaps(l: &mut [usize]) -> usize {
        let mut swaps: usize = 0;
        for i in (0..l.len()).rev() {
            for j in 0..i {
                if l[j] > l[j + 1] {
                    (l[j], l[j + 1]) = (l[j + 1], l[j]);
                    swaps += 1
                }
            }
        }
        swaps
    }

    // Generate list of objects in the algebra

    // The number of dimensions in the algebra
    // (e.g. use a 4D algebra to represent 3D geometry)
    let dimension = basis_vector_squares.len();

    let basis = {
        // We will represent a basis element in the algebra as a Vec<usize>
        // e.g. vec![2] = e_2, vec![1,2,3] = e_123, vec![] = 1

        // To generate all the basis elements, we will iterate through the k-vectors.
        // For each k, we right-multiply the set of basis vectors
        // onto the basis elements of grade k-1,
        // removing any results that are already
        // represented in the basis element set.

        // Start with the 0-vector, i.e. scalar, represented as vec![]
        let mut basis_km1_vectors = vec![vec![]];
        let mut basis = vec![vec![]];

        for _ in 1..=dimension {
            let mut basis_k_vectors = vec![];
            for b1 in basis_km1_vectors {
                for be2 in 0..dimension {
                    // Insert be2 into b1 in sorted order
                    match b1.binary_search(&be2) {
                        Ok(_) => {}
                        Err(pos) => {
                            let mut b = b1.clone();
                            b.insert(pos, be2);
                            if !basis_k_vectors.contains(&b) {
                                basis_k_vectors.push(b.clone());
                                basis.push(b);
                            }
                        }
                    }
                }
            }
            basis_km1_vectors = basis_k_vectors;
        }

        let basis_element_count = basis.len();

        // We now have a set of basis components in the vec `basis`.
        // Each one is a product of basis vectors, in sorted order.
        // This is a valid basis, but for improved ergonomics,
        // we will negate some of the basis elements to achieve two things:
        // 1. The ability to dualize an element in the geometry by reversing the list of coefficients
        //    e.g. in 4D algebra: 1e0 + 2e1 + 3e2 + 4e3 <-dual-> 4e012 + 3e031 + 2e023 + 1e132
        // 2. Maximize the number of basis elements that, when multiplied by their dual,
        //    equal I rather than -I. Note that we can make all of these products equal I in
        //    algebras with odd dimension, and we can always do at least half in even dimensions,
        //    but we are forced to choose between e.g. e123 or e132 in 4D
        //    i.e. choose between e0 * e123 = I and e132 * e0 = I
        //    (This code chooses the former--guaranteeing I over -I
        //    in the upper right quadrant of the multiplication table)

        for i in 0..basis_element_count {
            if i < basis_element_count / 2 {
                continue; // Do nothing to first half
            }

            let dual_i = basis_element_count - i - 1;

            let mut product: Vec<usize> = basis[dual_i]
                .iter()
                .cloned()
                .chain(basis[i].iter().cloned())
                .collect();
            let swaps = bubble_sort_count_swaps(product.as_mut());
            if swaps % 2 == 1 {
                // An odd number of swaps means that the product was equal to -I instead of I
                // and we must reverse the order of two of the elements
                let l = basis[i].len();
                (basis[i][l - 2], basis[i][l - 1]) = (basis[i][l - 1], basis[i][l - 2]);
            }
        }

        basis
    };

    fn coefficient_identifier(var: &str, basis: &[usize]) -> Ident {
        let basis_pretty: String = basis.iter().map(|e| format!("{}", e)).collect();
        Ident::new(&format!("{}{}", var, basis_pretty), Span::call_site())
    }

    ///////////////////////////////////////////////////////////////////////////

    // We will represent a multivector as an array of coefficients on the basis elements.
    // e.g. in 2D, there are 1 + 3 + 3 + 1 = 8 basis elements,
    // and a full multivector uses all of them: [1, 1, 1, 1, 1, 1, 1, 1]
    // An object such as a Bivector would only use a few of them: [0, 0, 0, 0, 1, 1, 1, 0]

    struct Object {
        name: Ident,
        select_components: DVector<bool>,
    }

    fn k_vector_identifier(k: usize, dimension: usize) -> Ident {
        Ident::new(
            &if k == dimension {
                "Pseudoscalar".to_owned()
            } else {
                match k {
                    0 => "Scalar".to_owned(),
                    1 => "Vector".to_owned(),
                    2 => "Bivector".to_owned(),
                    3 => "Trivector".to_owned(),
                    _ => format!("K{}Vector", k),
                }
            },
            Span::call_site(),
        )
    }

    let mut objects: Vec<Object> = Vec::new();

    // Add k-vectors to object set
    for k in 0..=dimension {
        let select_components =
            DVector::<bool>::from_iterator(basis.len(), basis.iter().map(|b| b.len() == k));
        let name = k_vector_identifier(k, dimension);
        objects.push(Object {
            name,
            select_components,
        });
    }

    // Add even & odd sub-alg"ebras to the object set
    for parity in 0..2 {
        let select_components = DVector::<bool>::from_iterator(
            basis.len(),
            basis.iter().map(|b| b.len() % 2 == parity),
        );
        let name = Ident::new(
            match parity {
                0 => "Even",
                1 => "Odd",
                _ => panic!("Expected parity to be 0 or 1"),
            },
            Span::call_site(),
        );
        objects.push(Object {
            name,
            select_components,
        });
    }

    // Generate struct & impl for each object
    let objects_code: TokenStream = objects
        .into_iter()
        .map(
            |Object {
                 name,
                 select_components,
             }| {
                let struct_members: TokenStream = basis
                    .iter()
                    .zip(select_components.iter())
                    .filter_map(|(b, is_selected)| {
                        is_selected.then(|| {
                            let identifier = coefficient_identifier("a", b);
                            quote! {
                                #identifier: T,
                            }
                        })
                    })
                    .collect();

                quote! {
                    // ===========================================================================
                    // #name
                    // ===========================================================================
                    struct #name<T> {
                        #struct_members
                    }
                }
            },
        )
        .collect();

    objects_code
}

#[proc_macro]
pub fn gen_algebra(input_tokens: proc_macro::TokenStream) -> proc_macro::TokenStream {
    // Parse input
    let input = parse_macro_input!(input_tokens as Input);

    gen_algebra2(input).into()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gen_algebra() {
        gen_algebra2(Input {
            basis_vector_squares: vec![1, 1, 1, 0],
        });
    }
}
