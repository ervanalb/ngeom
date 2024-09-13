extern crate proc_macro;
use core::cmp::Ordering;
use core::ops::{AddAssign, Mul};
use proc_macro2::Span;
use proc_macro2::{TokenStream, TokenTree};
use quote::{quote, ToTokens, TokenStreamExt};
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

    #[derive(Default, Clone)]
    struct SymbolicSumExpr(Vec<SymbolicProdExpr>);

    #[derive(PartialEq, Eq, Clone)]
    struct SymbolicProdExpr(isize, Vec<Symbol>);

    #[derive(PartialOrd, Ord, PartialEq, Eq, Clone)]
    struct Symbol(Ident);

    impl ToTokens for Symbol {
        fn to_tokens(&self, tokens: &mut TokenStream) {
            let Symbol(self_ident) = self;
            self_ident.to_tokens(tokens);
        }
    }

    impl PartialOrd for SymbolicProdExpr {
        fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
            Some(self.cmp(other))
        }
    }

    impl Ord for SymbolicProdExpr {
        fn cmp(&self, SymbolicProdExpr(other_coef, other_symbols): &Self) -> Ordering {
            let SymbolicProdExpr(self_coef, self_symbols) = self;
            self_symbols
                .cmp(&other_symbols)
                .then_with(|| self_coef.cmp(&other_coef))
        }
    }

    impl Mul<SymbolicProdExpr> for SymbolicProdExpr {
        type Output = SymbolicProdExpr;
        fn mul(mut self, SymbolicProdExpr(r_coef, mut r_symbols): Self) -> SymbolicProdExpr {
            let SymbolicProdExpr(l_coef, l_symbols) = &mut self;
            *l_coef *= r_coef;
            l_symbols.append(&mut r_symbols);
            self
        }
    }

    impl SymbolicProdExpr {
        fn simplify(mut self) -> Self {
            let SymbolicProdExpr(coef, symbols) = &mut self;
            // Sort expression
            if *coef == 0 {
                symbols.clear();
            } else {
                symbols.sort();
            }
            self
        }
    }

    impl ToTokens for SymbolicSumExpr {
        fn to_tokens(&self, tokens: &mut TokenStream) {
            let SymbolicSumExpr(terms) = self;
            if terms.len() == 0 {
                tokens.append_all(quote! { T::zero() });
            } else {
                for (count, prod_expr) in terms.iter().enumerate() {
                    let SymbolicProdExpr(coef, prod_terms) = prod_expr;
                    let coef = *coef;

                    if coef >= 0 {
                        if count != 0 {
                            tokens.append_all(quote! { + });
                        }
                    } else {
                        tokens.append_all(quote! { - });
                    }
                    let coef = coef.abs();

                    if prod_terms.len() == 0 {
                        // If there are no symbols in the product, then this is a scalar
                        if coef == 0 {
                            tokens.append_all(quote! { T::zero() });
                        } else if coef == 1 {
                            tokens.append_all(quote! {
                                T::one()
                            });
                        } else {
                            panic!("Scalar was not 0, -1 or 1");
                        }
                    } else {
                        // There are symbols in the product
                        if coef == 0 {
                            tokens.append_all(quote! { T::zero() * });
                        } else if coef == 1 {
                            // No token needed if coefficient is unity
                        } else if coef == 2 {
                            tokens.append_all(quote! { (T::one() + T::one()) * });
                        } else {
                            panic!("No representation for large coefficient {}", coef);
                        }
                        for (sym_count, sym) in prod_terms.iter().enumerate() {
                            if sym_count > 0 {
                                tokens.append_all(quote! { * });
                            }
                            sym.to_tokens(tokens);
                        }
                    }
                }
            }
        }
    }

    impl SymbolicSumExpr {
        fn simplify(self) -> Self {
            let SymbolicSumExpr(terms) = self;

            // Simplify all products
            let mut terms: Vec<_> = terms.into_iter().map(|prod| prod.simplify()).collect();

            // Sort expression by symbolic values
            terms.sort();

            // Combine adjacent terms whose symbolic parts are equal
            let mut new_expression = vec![];
            let mut prev_coef = 0;
            let mut prev_symbols = vec![];
            for SymbolicProdExpr(coef, symbols) in terms.into_iter() {
                if prev_symbols == symbols {
                    prev_coef += coef;
                } else {
                    new_expression.push(SymbolicProdExpr(prev_coef, prev_symbols));
                    prev_coef = coef;
                    prev_symbols = symbols;
                }
            }
            new_expression.push(SymbolicProdExpr(prev_coef, prev_symbols));

            let mut terms = new_expression;

            // Remove all products with coefficient = 0
            terms.retain(|SymbolicProdExpr(coef, _)| *coef != 0);

            SymbolicSumExpr(terms)
        }
    }

    impl AddAssign<SymbolicProdExpr> for SymbolicSumExpr {
        fn add_assign(&mut self, r_term: SymbolicProdExpr) {
            let SymbolicSumExpr(l_terms) = self;
            l_terms.push(r_term);
        }
    }

    impl AddAssign<SymbolicSumExpr> for SymbolicSumExpr {
        fn add_assign(&mut self, SymbolicSumExpr(mut r_terms): SymbolicSumExpr) {
            let SymbolicSumExpr(l_terms) = self;
            l_terms.append(&mut r_terms);
        }
    }

    impl Mul<SymbolicProdExpr> for SymbolicSumExpr {
        type Output = SymbolicSumExpr;
        fn mul(self, r: SymbolicProdExpr) -> SymbolicSumExpr {
            let SymbolicSumExpr(l) = self;
            SymbolicSumExpr(l.iter().map(|lp| lp.clone() * r.clone()).collect())
        }
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

        // We now have a set of basis components in the vec `basis`.
        // Each one is a product of basis vectors, in sorted order.

        // This is a valid basis, but for improved ergonomics,
        // we will negate some of the basis elements to minimize the number of negative signs in
        // the hodge star (dual) operator.
        // e.g. in 4D algebra: 1e0 + 2e1 + 3e2 <-dual-> 3e01 + 2e20 + 1e12
        // Note that we can make all of these products equal I in
        // algebras with odd dimension, and we can always do at least half in even dimensions,
        // but we are forced to choose between e.g. e123 or e132 in 4D
        // i.e. choose between e0 * e123 = I and e132 * e0 = I
        // (This code chooses the former--guaranteeing I over -I
        // in the upper right quadrant of the multiplication table)

        let basis_element_count = basis.len();

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

    let basis_element_count = basis.len();

    let right_complement_signs: Vec<_> = (0..basis_element_count)
        .map(|i| {
            let dual_i = basis_element_count - i - 1;

            // Compute the product of the basis element and its dual basis element
            // and figure out what the sign needs so that the product equals I
            // (rather than -I)
            let mut product: Vec<usize> = basis[i]
                .iter()
                .cloned()
                .chain(basis[dual_i].iter().cloned())
                .collect();
            let swaps = bubble_sort_count_swaps(product.as_mut());
            match swaps % 2 {
                0 => 1,
                1 => -1,
                _ => panic!("Expected parity to be 0 or 1"),
            }
        })
        .collect();

    // Returns the right_complement of a basis element as a pair of
    // (coef, complement_ix) where coef is -1 or 1
    fn right_complement(right_complement_signs: &Vec<isize>, i: usize) -> (isize, usize) {
        let complement_ix = right_complement_signs.len() - i - 1;
        (right_complement_signs[i], complement_ix)
    }

    // Returns the inverse of the right complement of a basis element
    // (i.e. the left complement)
    // as a pair of (coef, complement_ix) where coef is -1 or 1
    fn left_complement(right_complement_signs: &Vec<isize>, i: usize) -> (isize, usize) {
        let complement_ix = right_complement_signs.len() - i - 1;
        (right_complement_signs[complement_ix], complement_ix)
    }

    fn coefficient_ident(var: &str, basis: &[usize]) -> Ident {
        let basis_pretty: String = basis.iter().map(|e| format!("{}", e)).collect();
        Ident::new(&format!("{}{}", var, basis_pretty), Span::call_site())
    }

    // We will represent a multivector as an array of coefficients on the basis elements.
    // e.g. in 2D, there are 1 + 3 + 3 + 1 = 8 basis elements,
    // and a full multivector uses all of them: [1, 1, 1, 1, 1, 1, 1, 1]
    // An object such as a Bivector would only use a few of them: [0, 0, 0, 0, 1, 1, 1, 0]

    #[derive(PartialEq)]
    struct Object {
        name: Ident,
        select_components: Vec<bool>,
        is_scalar: bool,
        is_compound: bool,
    }

    impl Object {
        fn type_name(&self) -> TokenStream {
            if self.is_scalar {
                quote! { T }
            } else {
                let name = &self.name;
                quote! { #name < T > }
            }
        }
        fn type_name_colons(&self) -> TokenStream {
            if self.is_scalar {
                quote! { T }
            } else {
                let name = &self.name;
                quote! { #name :: < T > }
            }
        }
    }

    fn k_vector_identifier(k: usize, dimension: usize) -> Ident {
        Ident::new(
            &if k == dimension {
                "AntiScalar".to_owned()
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
        let select_components: Vec<_> = basis.iter().map(|b| b.len() == k).collect();
        let name = k_vector_identifier(k, dimension);
        let is_scalar = k == 0;
        objects.push(Object {
            name,
            select_components,
            is_scalar,
            is_compound: false,
        });
    }

    // Add scalar plus antiscalar to the object set
    let select_components: Vec<_> = basis
        .iter()
        .enumerate()
        .map(|(i, _)| i == 0 || i == basis.len() - 1)
        .collect();
    objects.push(Object {
        name: Ident::new("Magnitude", Span::call_site()),
        select_components,
        is_scalar: false,
        is_compound: true,
    });

    // Add even & odd sub-algebras to the object set
    for anti_parity in 0..2 {
        let anti_scalar_grade = basis[basis.len() - 1].len();
        let select_components: Vec<_> = basis
            .iter()
            .map(|b| (anti_scalar_grade - b.len()) % 2 == anti_parity)
            .collect();
        let name = Ident::new(
            match anti_parity {
                0 => "AntiEven",
                1 => "AntiOdd",
                _ => panic!("Expected parity to be 0 or 1"),
            },
            Span::call_site(),
        );
        objects.push(Object {
            name,
            select_components,
            is_scalar: false,
            is_compound: true,
        });
    }

    fn rewrite_identifiers<F: Fn(&Ident) -> Option<TokenStream>>(
        expressions_code: Vec<TokenStream>,
        replacement: F,
    ) -> Vec<TokenStream> {
        expressions_code
            .into_iter()
            .map(|expr| {
                expr.into_iter()
                    .flat_map(|tt| {
                        match &tt {
                            TokenTree::Ident(ident) => {
                                replacement(ident).unwrap_or(quote! { #ident })
                            }
                            other => quote! { #other },
                        }
                        .into_iter()
                    })
                    .collect()
            })
            .collect()
    }

    // Function to generate a simple unary operator on an object,
    // where entries can be rearranged and coefficients can be modified
    // (e.g. for implementing neg, right_complement, reverse)
    // The result will be a list of symbolic expressions,
    // one for each coefficient value in the resultant multivector.
    fn generate_symbolic_rearrangement<F: Fn(usize) -> (isize, usize)>(
        basis: &Vec<Vec<usize>>,
        select_components: &[bool],
        op: F,
    ) -> Vec<SymbolicSumExpr> {
        // Generate the sum
        let mut result: Vec<SymbolicSumExpr> = vec![Default::default(); select_components.len()];

        for (i, &is_selected) in select_components.iter().enumerate() {
            if is_selected {
                let (coef, result_basis_ix) = op(i);
                result[result_basis_ix] +=
                    SymbolicProdExpr(coef, vec![Symbol(coefficient_ident("a", &basis[i]))]);
            }
        }

        result.into_iter().map(|expr| expr.simplify()).collect()
    }

    fn generate_symbolic_norm<F: Fn(usize, usize) -> (isize, usize)>(
        basis: &Vec<Vec<usize>>,
        select_components: &[bool],
        product: F,
        sqrt: bool,
    ) -> Vec<SymbolicSumExpr> {
        // Generate the product
        let mut expressions: Vec<SymbolicSumExpr> =
            vec![Default::default(); select_components.len()];
        for i in select_components
            .iter()
            .enumerate()
            .filter_map(|(i, is_selected)| is_selected.then_some(i))
        {
            for j in select_components
                .iter()
                .enumerate()
                .filter_map(|(j, is_selected)| is_selected.then_some(j))
            {
                let (coef, ix_result) = product(i, j);

                expressions[ix_result] += SymbolicProdExpr(
                    coef,
                    vec![
                        Symbol(coefficient_ident("a", &basis[i])),
                        Symbol(coefficient_ident("a", &basis[j])),
                    ],
                );
            }
        }

        let expressions: Vec<_> = expressions
            .into_iter()
            .map(|expr| expr.simplify())
            .collect();

        if sqrt {
            // See if we can take the square root symbolically
            // Otherwise, return an empty expression (which will cause no code to be generated)
            {
                let is_scalar = expressions
                    .iter()
                    .enumerate()
                    .all(|(i, expr)| i == 0 || expr.0.len() == 0);
                let is_anti_scalar = expressions
                    .iter()
                    .enumerate()
                    .all(|(i, expr)| i == expressions.len() - 1 || expr.0.len() == 0);
                let expression = if is_scalar {
                    Some(expressions[0].clone())
                } else if is_anti_scalar {
                    Some(expressions[expressions.len() - 1].clone())
                } else {
                    None
                };
                if let Some(expression) = expression {
                    let SymbolicSumExpr(terms) = &expression;
                    if terms.len() == 1 {
                        let SymbolicProdExpr(coef, terms) = &terms[0];
                        if *coef == 1 && terms.len() == 2 && terms[0] == terms[1] {
                            let sqrt_expression =
                                SymbolicSumExpr(vec![SymbolicProdExpr(1, vec![terms[0].clone()])]);
                            let target_ix = if is_scalar {
                                0
                            } else if is_anti_scalar {
                                expressions.len() - 1
                            } else {
                                panic!("Took sqrt of something that wasn't a scalar or antiscalar");
                            };
                            Some(
                                (0..select_components.len())
                                    .map(|i| {
                                        if i == target_ix {
                                            sqrt_expression.clone()
                                        } else {
                                            Default::default()
                                        }
                                    })
                                    .collect(),
                            )
                        } else {
                            None // Expression is not a square
                        }
                    } else {
                        None // Multiple terms in the sum
                    }
                } else {
                    None // Squared norm is not a scalar or antiscalar
                }
            }
            .unwrap_or(vec![Default::default(); select_components.len()])
        } else {
            // Return the squared norm
            expressions
        }
    }

    fn gen_unary_operator(
        basis: &Vec<Vec<usize>>,
        objects: &[Object],
        op_trait: TokenStream,
        op_fn: Ident,
        obj: &Object,
        expressions: &[SymbolicSumExpr],
        alias: Option<(TokenStream, Ident)>,
    ) -> TokenStream {
        if obj.is_scalar {
            // Do not generate operations with the scalar being the LHS--
            // typically because these would violate rust's orphan rule
            // or result in conflicting trait implementations
            return quote! {};
        }

        // Figure out what the type of the output is
        let select_output_components: Vec<_> = expressions
            .iter()
            .map(|SymbolicSumExpr(e)| e.len() != 0)
            .collect();
        let output_object = objects.iter().find(|o| {
            o.select_components
                .iter()
                .zip(select_output_components.iter())
                .find(|&(obj_c, &out_c)| out_c && !obj_c)
                .is_none()
        });

        let Some(output_object) = output_object else {
            // No output object matches the result we got,
            // so don't generate any code
            return quote! {};
        };

        if output_object.is_scalar && expressions[0].0.len() == 0 {
            // This operation unconditionally returns 0,
            // so invoking it is probably a type error--
            // do not generate code for it
            return quote! {};
        }

        let output_type_name = &output_object.type_name_colons();
        let type_name = &obj.type_name();

        // Convert the expression to token trees
        let expressions: Vec<TokenStream> = expressions
            .into_iter()
            .map(|expr| quote! { #expr })
            .collect();

        // Find and replace temporary a### & b### identifier with self.a### & r.a###
        let expressions = rewrite_identifiers(expressions, |ident| {
            basis.iter().find_map(move |b| {
                if ident == &coefficient_ident("a", b) {
                    Some(if obj.is_scalar {
                        assert!(b.len() == 0);
                        quote! { self }
                    } else {
                        // All fields are in the format "a###"
                        // so we can re-use the temporary identifier
                        let new_ident = ident;
                        quote! { self.#new_ident }
                    })
                } else {
                    None
                }
            })
        });

        let return_expr = if output_object.is_scalar {
            let expr = &expressions[0];
            quote! { #expr }
        } else {
            let output_fields: TokenStream = output_object
                .select_components
                .iter()
                .zip(basis.iter().zip(expressions.iter()))
                .enumerate()
                .map(|(_i, (is_selected, (b, expr)))| {
                    if !is_selected {
                        return quote! {};
                    }
                    let field = coefficient_ident("a", b);
                    quote! { #field: #expr, }
                })
                .collect();
            quote! {
                #output_type_name {
                    #output_fields
                }
            }
        };

        let associated_output_type = quote! { type Output = #output_type_name; };

        let code = quote! {
            impl < T: Ring > #op_trait for #type_name {
                #associated_output_type

                fn #op_fn (self) -> #output_type_name {
                    #return_expr
                }
            }
        };

        let alias_code = if let Some((alias_trait, alias_fn)) = alias {
            quote! {
                impl < T: Ring > #alias_trait for #type_name {
                    #associated_output_type

                    fn #alias_fn (self) -> #output_type_name {
                        self.#op_fn()
                    }
                }
            }
        } else {
            quote! {}
        };

        quote! {
            #code
            #alias_code
        }
    }

    // Generate dot product multiplication table for the basis elements.
    // The dot product always produces a scalar.
    let dot_product_multiplication_table: Vec<Vec<isize>> = {
        let multiply_basis_vectors = |ei: usize, ej: usize| {
            // This would need to get more complicated for CGA
            // where the metric g is not a diagonal matrix
            if ei == ej {
                basis_vector_squares[ei]
            } else {
                0
            }
        };

        let multiply_basis_elements = |i: usize, j: usize| {
            // The scalar product of bivectors, trivectors, etc.
            // can be found using the Gram determinant.

            let bi = &basis[i];
            let bj = &basis[j];
            if bi.len() != bj.len() {
                return 0;
            }
            if bi.len() == 0 {
                return 1; // 1 • 1 = 1
            }

            let gram_matrix: Vec<Vec<isize>> = bi
                .iter()
                .map(|&ei| {
                    bj.iter()
                        .map(|&ej| multiply_basis_vectors(ei, ej))
                        .collect()
                })
                .collect();

            fn determinant(m: &Vec<Vec<isize>>) -> isize {
                if m.len() == 1 {
                    m[0][0]
                } else {
                    let n = m.len();
                    (0..n)
                        .map(move |j| {
                            let i = 0;
                            let sign = match (i + j) % 2 {
                                0 => 1,
                                1 => -1,
                                _ => panic!("Expected parity to be 0 or 1"),
                            };

                            let minor: Vec<Vec<_>> = (0..n)
                                .flat_map(|i2| {
                                    if i2 == i {
                                        None
                                    } else {
                                        Some(
                                            (0..n)
                                                .flat_map(|j2| {
                                                    if j2 == j {
                                                        None
                                                    } else {
                                                        Some(m[i2][j2])
                                                    }
                                                })
                                                .collect(),
                                        )
                                    }
                                })
                                .collect();

                            sign * m[i][j] * determinant(&minor)
                        })
                        .sum()
                }
            }
            determinant(&gram_matrix)
        };

        (0..basis_element_count)
            .map(|i| {
                (0..basis_element_count)
                    .map(move |j| multiply_basis_elements(i, j))
                    .collect()
            })
            .collect()
    };

    // Generate geometric product multiplication table for the basis elements
    // Each entry is a tuple of the (coefficient, basis_index)
    // e.g. (1, 0) means the multiplication result is 1 * scalar = 1
    let geometric_product_multiplication_table: Vec<Vec<(isize, usize)>> = {
        let multiply_basis_elements = |i: usize, j: usize| {
            let mut product: Vec<_> = basis[i]
                .iter()
                .cloned()
                .chain(basis[j].iter().cloned())
                .collect();
            let swaps = bubble_sort_count_swaps(product.as_mut());
            let mut coef = match swaps % 2 {
                0 => 1,
                1 => -1,
                _ => panic!("Expected parity to be 0 or 1"),
            };

            // Remove repeated elements in the product
            let mut new_product = vec![];
            let mut prev_e = None;
            for e in product.into_iter() {
                if Some(e) == prev_e {
                    coef *= basis_vector_squares[e];
                    prev_e = None;
                } else {
                    if let Some(prev_e) = prev_e {
                        new_product.push(prev_e);
                    }
                    prev_e = Some(e);
                }
            }
            if let Some(prev_e) = prev_e {
                new_product.push(prev_e);
            }

            // Figure out which basis element this corresponds to
            basis
                .iter()
                .enumerate()
                .find_map(|(i, b)| {
                    let mut b_sorted = b.clone();
                    let swaps = bubble_sort_count_swaps(b_sorted.as_mut());
                    (new_product == b_sorted).then(|| {
                        let coef = coef
                            * match swaps % 2 {
                                0 => 1,
                                1 => -1,
                                _ => panic!("Expected parity to be 0 or 1"),
                            };
                        (coef, i)
                    })
                })
                .expect("Product of basis elements not found in basis set")
        };

        (0..basis_element_count)
            .map(|i| {
                (0..basis_element_count)
                    .map(move |j| multiply_basis_elements(i, j))
                    .collect()
            })
            .collect()
    };

    let geometric_anti_product_multiplication_table: Vec<Vec<(isize, usize)>> = {
        let multiply_basis_elements = |i: usize, j: usize| {
            let (coef_i, i) = right_complement(&right_complement_signs, i);
            let (coef_j, j) = right_complement(&right_complement_signs, j);
            let (coef, ix) = geometric_product_multiplication_table[i][j];
            let (coef_comp, ix) = left_complement(&right_complement_signs, ix);
            (coef_i * coef_j * coef * coef_comp, ix)
        };

        (0..basis_element_count)
            .map(|i| {
                (0..basis_element_count)
                    .map(move |j| multiply_basis_elements(i, j))
                    .collect()
            })
            .collect()
    };

    // Function to generate a sum of two objects, e.g. for overloading + or -
    // The result will be a list of symbolic expressions,
    // one for each coefficient value in the resultant multivector.
    fn generate_symbolic_sum(
        basis: &Vec<Vec<usize>>,
        select_components_a: &[bool],
        select_components_b: &[bool],
        coef_a: isize,
        coef_b: isize,
    ) -> Vec<SymbolicSumExpr> {
        // Generate the sum
        let mut result: Vec<SymbolicSumExpr> = vec![Default::default(); select_components_a.len()];

        for (i, (&a_is_selected, &b_is_selected)) in select_components_a
            .iter()
            .zip(select_components_b.iter())
            .enumerate()
        {
            if a_is_selected {
                result[i] +=
                    SymbolicProdExpr(coef_a, vec![Symbol(coefficient_ident("a", &basis[i]))]);
            }
            if b_is_selected {
                result[i] +=
                    SymbolicProdExpr(coef_b, vec![Symbol(coefficient_ident("b", &basis[i]))]);
            }
        }

        result.into_iter().map(|expr| expr.simplify()).collect()
    }

    // Function to generate a product of two objects, e.g. geometric product, wedge product, etc.
    // The result will be a list of symbolic expressions,
    // one for each coefficient value in the resultant multivector.
    fn generate_symbolic_product<F: Fn(usize, usize) -> (isize, usize)>(
        basis: &Vec<Vec<usize>>,
        select_components_a: &[bool],
        select_components_b: &[bool],
        product: F,
    ) -> Vec<SymbolicSumExpr> {
        // Generate the product
        let mut result: Vec<SymbolicSumExpr> = vec![Default::default(); select_components_a.len()];
        for i in select_components_a
            .iter()
            .enumerate()
            .filter_map(|(i, is_selected)| is_selected.then_some(i))
        {
            for j in select_components_b
                .iter()
                .enumerate()
                .filter_map(|(j, is_selected)| is_selected.then_some(j))
            {
                let (coef, result_basis_ix) = product(i, j);

                result[result_basis_ix] += SymbolicProdExpr(
                    coef,
                    vec![
                        Symbol(coefficient_ident("a", &basis[i])),
                        Symbol(coefficient_ident("b", &basis[j])),
                    ],
                );
            }
        }

        result.into_iter().map(|expr| expr.simplify()).collect()
    }

    // Function to generate a double product of two objects, e.g. sandwich product, project,
    // etc.
    // The result will be a list of symbolic expressions,
    // one for each coefficient value in the resultant multivector.
    // The resulting code will implement the product in the following order:
    // (B PRODUCT1 A) PRODUCT2 B
    fn generate_symbolic_double_product<
        F1: Fn(usize, usize) -> (isize, usize),
        F2: Fn(usize, usize) -> (isize, usize),
    >(
        basis: &Vec<Vec<usize>>,
        select_components_a: &[bool],
        select_components_b: &[bool],
        product_1: F1,
        product_2: F2,
    ) -> Vec<SymbolicSumExpr> {
        // Generate the first intermediate product B PRODUCT1 A
        // where i maps to components of B, and j maps to components of A
        let mut intermediate_result: Vec<SymbolicSumExpr> =
            vec![Default::default(); select_components_a.len()];
        for i in select_components_b
            .iter()
            .enumerate()
            .filter_map(|(i, is_selected)| is_selected.then_some(i))
        {
            for j in select_components_a
                .iter()
                .enumerate()
                .filter_map(|(j, is_selected)| is_selected.then_some(j))
            {
                let (coef, result_basis_ix) = product_1(i, j);
                intermediate_result[result_basis_ix] += SymbolicProdExpr(
                    coef,
                    vec![
                        Symbol(coefficient_ident("b", &basis[i])),
                        Symbol(coefficient_ident("a", &basis[j])),
                    ],
                );
            }
        }
        let intermediate_result: Vec<_> = intermediate_result
            .into_iter()
            .map(|expr| expr.simplify())
            .collect();

        // Generate the final product (B PRODUCT1 A) PRODUCT2 B
        // where i maps to components of the intermediate result B PRODUCT1 A
        // and j maps to components of B.
        let mut result: Vec<SymbolicSumExpr> = vec![Default::default(); select_components_a.len()];
        for (i, intermediate_term) in intermediate_result.iter().enumerate() {
            for j in select_components_b
                .iter()
                .enumerate()
                .filter_map(|(j, is_selected)| is_selected.then_some(j))
            {
                let (coef, result_basis_ix) = product_2(i, j);
                let new_term =
                    SymbolicProdExpr(coef, vec![Symbol(coefficient_ident("b", &basis[j]))]);
                let result_term = intermediate_term.clone() * new_term;
                result[result_basis_ix] += result_term;
            }
        }

        result.into_iter().map(|expr| expr.simplify()).collect()
    }

    fn gen_binary_operator<F: Fn(&[bool], &[bool]) -> Vec<SymbolicSumExpr>>(
        basis: &Vec<Vec<usize>>,
        objects: &[Object],
        op_trait: TokenStream,
        op_fn: Ident,
        lhs_obj: &Object,
        op: F,
        implicit_promotion_to_compound: bool,
        alias: Option<(TokenStream, Ident)>,
    ) -> TokenStream {
        objects
            .iter()
            .map(|rhs_obj| {
                if lhs_obj.is_scalar {
                    // Do not generate operations with the scalar being the LHS--
                    // these violate rust's orphan rule
                    // or result in conflicting trait implementations.
                    // Technically we could allow this for custom ops such as T::cross(r: Vector<T>)
                    // but we won't for the sake of consistency.
                    // Use Vector<T>::cross(r: T) instead.
                    return quote! {};
                }

                let expressions = op(&lhs_obj.select_components, &rhs_obj.select_components);

                let select_output_components: Vec<_> = expressions
                    .iter()
                    .map(|SymbolicSumExpr(e)| e.len() != 0)
                    .collect();
                let output_object = objects.iter().find(|o| {
                    o.select_components
                        .iter()
                        .zip(select_output_components.iter())
                        .find(|&(obj_c, &out_c)| out_c && !obj_c)
                        .is_none()
                });

                let Some(output_object) = output_object else {
                    // No output object matches the result we got,
                    // so don't generate any code
                    return quote! {};
                };

                let rhs_type_name = &rhs_obj.type_name();
                let lhs_type_name = &lhs_obj.type_name();
                let output_type_name = &output_object.type_name_colons();

                if !implicit_promotion_to_compound
                    && (output_object.is_compound
                        && !(lhs_obj.is_scalar
                            || rhs_obj.is_scalar
                            || (lhs_obj.is_compound && rhs_obj.is_compound)))
                {
                    // Do not create compound objects unintentionally. Only allow:
                    // 1. Products of compound objects and scalars
                    // 2. Products of compound objects and other compound objects
                    return quote! {};
                }

                if output_object.is_scalar && expressions[0].0.len() == 0 {
                    // This operation unconditionally returns 0,
                    // so invoking it is probably a type error--
                    // do not generate code for it.

                    return quote! {};
                }

                // Convert the expression to token trees
                let expressions: Vec<TokenStream> = expressions
                    .into_iter()
                    .map(|expr| quote! { #expr })
                    .collect();

                // Find and replace temporary a### & b### identifier with self.a### & r.a###
                let expressions = rewrite_identifiers(expressions, |ident| {
                    basis.iter().find_map(move |b| {
                        if ident == &coefficient_ident("a", b) {
                            Some(if lhs_obj.is_scalar {
                                assert!(b.len() == 0);
                                quote! { self }
                            } else {
                                // All fields are in the format "a###"
                                // so we can re-use the temporary identifier
                                let new_ident = ident;
                                quote! { self.#new_ident }
                            })
                        } else if ident == &coefficient_ident("b", b) {
                            Some(if rhs_obj.is_scalar {
                                assert!(b.len() == 0);
                                quote! { r }
                            } else {
                                // All fields are in the format "a###"
                                // so we need to convert b### to a###
                                let new_ident = coefficient_ident("a", b);
                                quote! { r.#new_ident }
                            })
                        } else {
                            None
                        }
                    })
                });

                let return_expr = if output_object.is_scalar {
                    let expr = &expressions[0];
                    quote! { #expr }
                } else {
                    let output_fields: TokenStream = output_object
                        .select_components
                        .iter()
                        .zip(basis.iter().zip(expressions.iter()))
                        .enumerate()
                        .map(|(_i, (is_selected, (b, expr)))| {
                            if !is_selected {
                                return quote! {};
                            }
                            let field = coefficient_ident("a", b);
                            quote! { #field: #expr, }
                        })
                        .collect();
                    quote! {
                        #output_type_name {
                            #output_fields
                        }
                    }
                };

                let associated_output_type = quote! { type Output = #output_type_name; };

                let code = quote! {
                    impl < T: Ring > #op_trait < #rhs_type_name >  for #lhs_type_name {
                        #associated_output_type

                        fn #op_fn (self, r: #rhs_type_name) -> #output_type_name {
                            #return_expr
                        }
                    }
                };

                let alias_code = if let Some((alias_trait, alias_fn)) = &alias {
                    quote! {
                        impl < T: Ring > #alias_trait < #rhs_type_name > for #lhs_type_name {
                            #associated_output_type

                            fn #alias_fn (self, r: #rhs_type_name) -> #output_type_name {
                                self.#op_fn(r)
                            }
                        }
                    }
                } else {
                    quote! {}
                };

                quote! {
                    #code
                    #alias_code
                }
            })
            .collect()
    }

    // Generate struct & impl for each object
    let struct_code: TokenStream = objects
        .iter()
        .map(|obj| {
            // The struct definition
            let struct_definition = match obj.is_scalar {
                true => quote! {}, // Use T for scalar type--don't wrap it in a struct
                false => {
                    let struct_members: TokenStream = basis
                        .iter()
                        .zip(obj.select_components.iter())
                        .filter_map(|(b, is_selected)| {
                            is_selected.then(|| {
                                let identifier = coefficient_ident("a", b);
                                quote! {
                                    pub #identifier: T,
                                }
                            })
                        })
                        .collect();
                    let name = &obj.name;
                    quote! {
                        pub struct #name<T> {
                            #struct_members
                        }
                    }
                }
            };

            quote! {
                #struct_definition
            }
        })
        .collect();

    let anti_scalar_ops: TokenStream = {
        let field = coefficient_ident("a", &basis[basis.len() - 1]);

        quote! {
            impl<T: Abs> AntiAbs for AntiScalar<T> {
                type Output = AntiScalar<<T as Abs>::Output>;
                fn anti_abs(self) -> Self::Output {
                    AntiScalar {
                        #field: self.#field.abs(),
                    }
                }
            }

            impl<T: Recip> AntiRecip for AntiScalar<T> {
                type Output = AntiScalar<<T as Recip>::Output>;
                fn anti_recip(self) -> Self::Output {
                    AntiScalar {
                        #field: self.#field.recip(),
                    }
                }
            }

            impl<T: Sqrt> AntiSqrt for AntiScalar<T> {
                fn anti_sqrt(self) -> AntiScalar<T> {
                    AntiScalar {
                        #field: self.#field.sqrt(),
                    }
                }
            }

            impl<T: Trig> AntiTrig for AntiScalar<T> {
                type Output = AntiScalar<<T as Trig>::Output>;

                fn anti_cos(self) -> Self::Output {
                    AntiScalar {
                        #field: self.#field.cos(),
                    }
                }
                fn anti_sin(self) -> Self::Output {
                    AntiScalar {
                        #field: self.#field.sin(),
                    }
                }
                fn anti_sinc(self) -> Self::Output {
                    AntiScalar {
                        #field: self.#field.sinc(),
                    }
                }
            }

            impl<T: core::cmp::PartialOrd> PartialOrd for AntiScalar<T> {
                fn partial_cmp(&self, other: &Self) -> Option<core::cmp::Ordering> {
                    self.#field.partial_cmp(&other.#field)
                }
            }

            impl<T: core::cmp::Ord> Ord for AntiScalar<T> {
                fn cmp(&self, other: &Self) -> core::cmp::Ordering {
                    self.#field.cmp(&other.#field)
                }
            }
        }
    };

    let impl_code: TokenStream = objects
        .iter()
        .map(|obj| {
            // Derive clone
            let clone_code = {
                match obj.is_scalar {
                    true => quote! {}, // Use T for scalar type--don't wrap it in a struct
                    false => {
                        let type_name = &obj.type_name();
                        let output_type_name = &obj.type_name_colons();

                        let cloned_fields: TokenStream = basis
                            .iter()
                            .zip(obj.select_components.iter())
                            .filter_map(|(b, is_selected)| {
                                is_selected.then(|| {
                                    let identifier = coefficient_ident("a", b);
                                    quote! {
                                        #identifier: self.#identifier.clone(),
                                    }
                                })
                            })
                            .collect();
                        quote! {
                            impl <T: core::clone::Clone> core::clone::Clone for #type_name {
                                fn clone(&self) -> Self {
                                    #output_type_name {
                                        #cloned_fields
                                    }
                                }
                            }
                        }
                    }
                }
            };
            // Derive copy
            let copy_code = {
                match obj.is_scalar {
                    true => quote! {}, // Use T for scalar type--don't wrap it in a struct
                    false => {
                        let type_name = &obj.type_name();
                        quote! {
                            impl <T: core::marker::Copy> core::marker::Copy for #type_name { }
                        }
                    }
                }
            };
            
            // Derive debug
            let debug_code = {
                match obj.is_scalar {
                    true => quote! {}, // Use T for scalar type--don't wrap it in a struct
                    false => {
                        let fmt_expr_components: String = basis
                            .iter()
                            .zip(obj.select_components.iter())
                            .filter_map(|(b, is_selected)| {
                                is_selected.then(||
                                    format!("{}: {{:?}}, ", coefficient_ident("a", b))
                                )
                            })
                            .collect();
                        let fmt_vars: TokenStream = basis
                            .iter()
                            .zip(obj.select_components.iter())
                            .filter_map(|(b, is_selected)| {
                                is_selected.then(|| {
                                    let identifier = coefficient_ident("a", b);
                                    quote! {
                                        self.#identifier,
                                    }
                                })
                            })
                            .collect();
                        let fmt_expr = format!("{}: {{{{{}}}}}", obj.name, fmt_expr_components);
                        let type_name = &obj.type_name_colons();
                        quote! {
                            impl <T: core::fmt::Debug> core::fmt::Debug for #type_name {
                                fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
                                    write!(f, #fmt_expr, #fmt_vars)
                                }
                            }
                        }
                    }
                }
            };

            // Derive display
            let display_code = {
                match obj.is_scalar {
                    true => quote! {}, // Use T for scalar type--don't wrap it in a struct
                    false => {
                        let fmt_expr: String = basis
                            .iter()
                            .zip(obj.select_components.iter())
                            .filter_map(|(b, is_selected)| {
                                is_selected.then(||
                                    if b.len() == 0 {
                                        "({}) + ".to_string()
                                    } else {
                                        format!("({{}}) * {} + ", coefficient_ident("e", b))
                                    }
                                )
                            })
                            .collect();
                        let fmt_vars: TokenStream = basis
                            .iter()
                            .zip(obj.select_components.iter())
                            .filter_map(|(b, is_selected)| {
                                is_selected.then(|| {
                                    let identifier = coefficient_ident("a", b);
                                    quote! {
                                        self.#identifier,
                                    }
                                })
                            })
                            .collect();
                        let type_name = &obj.type_name_colons();
                        quote! {
                            impl <T: core::fmt::Display> core::fmt::Display for #type_name {
                                fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
                                    write!(f, #fmt_expr, #fmt_vars)
                                }
                            }
                        }
                    }
                }
            };

            // Derive PartialEq
            let partial_eq_code = {
                match obj.is_scalar {
                    true => quote! {}, // Use T for scalar type--don't wrap it in a struct
                    false => {
                        let type_name = &obj.type_name();

                        let eq_fields: TokenStream = basis
                            .iter()
                            .zip(obj.select_components.iter())
                            .filter_map(|(b, is_selected)| is_selected.then_some(b))
                            .enumerate()
                            .map(|(i, b)| {
                                let and = match i {
                                    0 => quote! {},
                                    _ => quote! { && },
                                };
                                let identifier = coefficient_ident("a", b);
                                quote! {
                                    #and self.#identifier == other.#identifier
                                }
                            })
                            .collect();
                        quote! {
                            impl <T: core::cmp::PartialEq> core::cmp::PartialEq for #type_name {
                                fn eq(&self, other: &Self) -> bool {
                                    #eq_fields
                                }
                            }
                        }
                    }
                }
            };

            // Derive Eq
            let eq_code = {
                match obj.is_scalar {
                    true => quote! {}, // Use T for scalar type--don't wrap it in a struct
                    false => {
                        let type_name = &obj.type_name();
                        quote! {
                            impl <T: core::cmp::Eq> core::cmp::Eq for #type_name { }
                        }
                    }
                }
            };

            // Overload unary -
            let op_trait = quote! { core::ops::Neg };
            let op_fn = Ident::new("neg", Span::call_site());
            let neg_code = gen_unary_operator(
                &basis,
                &objects,
                op_trait,
                op_fn,
                &obj,
                &generate_symbolic_rearrangement(&basis, &obj.select_components, |i: usize| (-1, i)),
                None,
            );

            // Add a method A.reverse()
            let op_trait = quote! { Reverse };
            let op_fn = Ident::new("reverse", Span::call_site());
            let reverse_f = |i: usize| {
                    let coef = match (basis[i].len() / 2) % 2 {
                        0 => 1,
                        1 => -1,
                        _ => panic!("Expected parity to be 0 or 1"),
                    };
                    (coef, i)
                };
            let reverse_expressions = generate_symbolic_rearrangement(&basis, &obj.select_components, reverse_f);
            let reverse_code = gen_unary_operator(&basis, &objects, op_trait, op_fn, &obj, &reverse_expressions, None);

            // Add a method A.anti_reverse()
            let op_trait = quote! { AntiReverse };
            let op_fn = Ident::new("anti_reverse", Span::call_site());
            let anti_reverse_f = |i: usize| {
                    let (coef_i, i) = right_complement(&right_complement_signs, i);
                    let (coef_rev, ix) = reverse_f(i);
                    let (coef_comp, ix) = left_complement(&right_complement_signs, ix);
                    (coef_i * coef_rev * coef_comp, ix)
                };
            let anti_reverse_expressions = generate_symbolic_rearrangement(&basis, &obj.select_components, anti_reverse_f);
            let anti_reverse_code = gen_unary_operator(&basis, &objects, op_trait, op_fn, &obj, &anti_reverse_expressions, None);

            // Add a method A.right_complement()
            let op_trait = quote! { RightComplement };
            let op_fn = Ident::new("right_complement", Span::call_site());
            let right_complement_expressions = generate_symbolic_rearrangement(&basis, &obj.select_components, |i: usize| right_complement(&right_complement_signs, i));
            let right_complement_code = gen_unary_operator(&basis, &objects, op_trait, op_fn, &obj, &right_complement_expressions, None);

            // Add a method A.left_complement()
            let op_trait = quote! { LeftComplement };
            let op_fn = Ident::new("left_complement", Span::call_site());
            let left_complement_expressions = generate_symbolic_rearrangement(&basis, &obj.select_components, |i: usize| left_complement(&right_complement_signs, i));
            let left_complement_code = gen_unary_operator(&basis, &objects, op_trait, op_fn, &obj, &left_complement_expressions, None);

            // Add a method A.bulk()
            let op_trait = quote! { Bulk };
            let op_fn = Ident::new("bulk", Span::call_site());
            let bulk_f = |i: usize| {
                // This would need to be changed for CGA where the dot product multiplication table
                // is not diagonal
                (dot_product_multiplication_table[i][i], i)
            };
            let bulk_expressions = generate_symbolic_rearrangement(&basis, &obj.select_components, bulk_f);
            let bulk_code = gen_unary_operator(&basis, &objects, op_trait, op_fn, &obj, &bulk_expressions, None);

            // Add a method A.weight()
            let op_trait = quote! { Weight };
            let op_fn = Ident::new("weight", Span::call_site());
            let weight_f = |i: usize| {
                let (coef_i, i) = right_complement(&right_complement_signs, i);
                let (coef_metric, ix) = bulk_f(i);
                let (coef_complement, ix) = left_complement(&right_complement_signs, ix);
                (coef_i * coef_metric * coef_complement, ix)
            };
            let weight_expressions = generate_symbolic_rearrangement(&basis, &obj.select_components, weight_f);
            let weight_code = gen_unary_operator(&basis, &objects, op_trait, op_fn, &obj, &weight_expressions, None);

            // Add a method A.bulk_dual() which computes A★
            let op_trait = quote! { BulkDual };
            let op_fn = Ident::new("bulk_dual", Span::call_site());
            let bulk_dual_f = |i: usize| {
                let (coef_bulk, i) = bulk_f(i);
                let (coef_comp, i) = right_complement(&right_complement_signs, i);
                (coef_bulk * coef_comp, i)
            };
            let bulk_dual_expressions = generate_symbolic_rearrangement(&basis, &obj.select_components, bulk_dual_f);
            let bulk_dual_code =
                gen_unary_operator(&basis, &objects, op_trait, op_fn, &obj, &bulk_dual_expressions, None);

            // Add a method A.weight_dual() which computes A☆
            let op_trait = quote! { WeightDual };
            let op_fn = Ident::new("weight_dual", Span::call_site());
            let alias_trait = quote! { Normal };
            let alias_fn = Ident::new("normal", Span::call_site());
            let weight_dual_f = |i: usize| {
                let (coef_bulk, i) = weight_f(i);
                let (coef_comp, i) = right_complement(&right_complement_signs, i);
                (coef_bulk * coef_comp, i)
            };
            let weight_dual_expressions = generate_symbolic_rearrangement(&basis, &obj.select_components, weight_dual_f);
            let weight_dual_code =
                gen_unary_operator(&basis, &objects, op_trait, op_fn, &obj, &weight_dual_expressions, Some((alias_trait, alias_fn)));

            // Add a method A.bulk_norm_squared()
            let op_trait = quote! { BulkNormSquared };
            let op_fn = Ident::new("bulk_norm_squared", Span::call_site());
            // Squared norm uses the product A • A
            let bulk_norm_product_f = |i: usize, j: usize| (dot_product_multiplication_table[i][j], 0);

            let bulk_norm_squared_expressions = generate_symbolic_norm(&basis, &obj.select_components, bulk_norm_product_f, false);
            let bulk_norm_squared_code =
                gen_unary_operator(&basis, &objects, op_trait, op_fn, &obj, &bulk_norm_squared_expressions, None);

            let bulk_norm_code = if !bulk_norm_squared_code.is_empty() {
                // Add a method A.bulk_norm() if it is possible to symbolically take the sqrt
                let op_trait = quote! { BulkNorm };
                let op_fn = Ident::new("bulk_norm", Span::call_site());
                let bulk_norm_expressions = generate_symbolic_norm(&basis, &obj.select_components, bulk_norm_product_f, true);
                let bulk_norm_code = gen_unary_operator(&basis, &objects, op_trait, op_fn, &obj, &bulk_norm_expressions, None);
                if !bulk_norm_code.is_empty() {
                    // We found a way to symbolically take the sqrt
                    bulk_norm_code
                } else {
                    // Take the square root of the norm numerically
                    let type_name = obj.type_name();
                    quote! {
                        impl < T: Ring + Sqrt > BulkNorm for #type_name
                            where <Self as BulkNormSquared>::Output: Sqrt {
                            type Output = <Self as BulkNormSquared>::Output;

                            fn bulk_norm (self) -> Self::Output {
                                self.bulk_norm_squared().sqrt()
                            }
                        }
                    }
                }
            } else {
                quote! {}  // There is no squared bulk norm, so return no impls for norm
            };

            let op_trait = quote! { WeightNormSquared };
            let op_fn = Ident::new("weight_norm_squared", Span::call_site());
            let weight_norm_product_f = |i: usize, j: usize| {
                let (coef_i, i) = right_complement(&right_complement_signs, i);
                let (coef_j, j) = right_complement(&right_complement_signs, j);
                let (coef_prod, ix) = bulk_norm_product_f(i, j);
                let (coef_complement, ix) = left_complement(&right_complement_signs, ix);
                (coef_i * coef_j * coef_prod * coef_complement, ix)
            };
            let weight_norm_squared_expressions = generate_symbolic_norm(&basis, &obj.select_components, weight_norm_product_f, false);

            let weight_norm_squared_code =
                gen_unary_operator(&basis, &objects, op_trait, op_fn, &obj, &weight_norm_squared_expressions, None);

            let weight_norm_code = if !weight_norm_squared_code.is_empty() {
                // Add a method A.weight_norm() if it is possible to symbolically take the sqrt
                let op_trait = quote! { WeightNorm };
                let op_fn = Ident::new("weight_norm", Span::call_site());
                let weight_norm_expressions = generate_symbolic_norm(&basis, &obj.select_components, weight_norm_product_f, true);
                let weight_norm_code =
                    gen_unary_operator(&basis, &objects, op_trait, op_fn, &obj, &weight_norm_expressions, None);
                if !weight_norm_code.is_empty() {
                    // We found a way to symbolically take the sqrt
                    weight_norm_code
                } else {
                    // Take the square root of the norm numerically
                    let type_name = obj.type_name();
                    quote! {
                        impl < T: Ring + Sqrt > WeightNorm for #type_name
                            where <Self as WeightNormSquared>::Output: AntiSqrt {
                            type Output = <Self as WeightNormSquared>::Output;

                            fn weight_norm (self) -> Self::Output {
                                self.weight_norm_squared().anti_sqrt()
                            }
                        }
                    }
                }
            } else {
                quote! {}  // There is no squared norm, so return no impls for norm
            };

            /*
            // Add a method A.geometric_norm_squared()
            let op_trait = quote! { GeometricNormSquared };
            let op_fn = Ident::new("geometric_norm_squared", Span::call_site());
            let geometric_norm_product_f =|i: usize, j: usize| {
                // Compute the product A * ~A

                let reverse_coef = match (basis[j].len() / 2) % 2 {
                    0 => 1,
                    1 => -1,
                    _ => panic!("Expected parity to be 0 or 1"),
                };

                let (coef, ix_result) = multiplication_table[i][j];

                (coef * reverse_coef, ix_result)
            };
            let geometric_norm_squared_expressions = generate_symbolic_norm(&basis, &obj.select_components, geometric_norm_product_f, false);
            let geometric_norm_squared_code =
                gen_unary_operator(&basis, &objects, op_trait, op_fn, &obj, &geometric_norm_squared_expressions
                    , false);

            let geometric_norm_code = if !geometric_norm_squared_code.is_empty() {
                // Add a method A.geometric_norm() if it is possible to symbolically take the sqrt
                let op_trait = quote! { GeometricNorm };
                let op_fn = Ident::new("geometric_norm", Span::call_site());
                let geometric_norm_code =
                    gen_unary_operator(&basis, &objects, op_trait, op_fn, &obj, &generate_symbolic_norm(&basis, &obj.select_components, geometric_norm_product_f, true)
                        , false);
                if !geometric_norm_code.is_empty() {
                    // We found a way to symbolically take the sqrt
                    geometric_norm_code
                } else {
                    // Take the square root of the norm numerically
                    let type_name = obj.type_name();
                    quote! {
                        impl < T: Ring > GeometricNorm for #type_name
                            where <Self as GeometricNormSquared>::Output: Sqrt {
                            type Output = <Self as GeometricNormSquared>::Output;

                            fn geometric_norm (self) -> Self::Output {
                                self.geometric_norm_squared().sqrt()
                            }
                        }
                    }
                }
            } else {
                quote! {}  // There is no squared norm, so return no impls for norm
            };

            let op_trait = quote! { GeometricInfNormSquared };
            let op_fn = Ident::new("geometric_inf_norm_squared", Span::call_site());
            let geometric_inf_norm_product_f = |i: usize, j: usize| {
                let (coef_i, i) = right_complement(&right_complement_signs, i);
                let (coef_j, j) = right_complement(&right_complement_signs, j);
                let (coef_prod, ix) = geometric_norm_product_f(i, j);
                (coef_i * coef_j * coef_prod, ix)
            };
            let geometric_inf_norm_squared_expressions = generate_symbolic_norm(&basis, &obj.select_components, geometric_inf_norm_product_f, false);
            let geometric_inf_norm_squared_code =
                gen_unary_operator(&basis, &objects, op_trait, op_fn, &obj, &geometric_inf_norm_squared_expressions
                    , false);

            let geometric_inf_norm_code = if !geometric_inf_norm_squared_code.is_empty() {
                // Add a method A.geometric_inf_norm() if it is possible to symbolically take the sqrt
                let op_trait = quote! { GeometricInfNorm };
                let op_fn = Ident::new("geometric_inf_norm", Span::call_site());
                let geometric_inf_norm_code =
                    gen_unary_operator(&basis, &objects, op_trait, op_fn, &obj, &generate_symbolic_norm(&basis, &obj.select_components, geometric_inf_norm_product_f, true)
                        , false);
                if !geometric_inf_norm_code.is_empty() {
                    // We found a way to symbolically take the sqrt
                    geometric_inf_norm_code
                } else {
                    // Take the square root of the norm numerically
                    let type_name = obj.type_name();
                    quote! {
                        impl < T: Ring > GeometricInfNorm for #type_name
                            where <Self as GeometricInfNormSquared>::Output: Sqrt {
                            type Output = <Self as GeometricInfNormSquared>::Output;

                            fn geometric_inf_norm (self) -> Self::Output {
                                self.geometric_inf_norm_squared().sqrt()
                            }
                        }
                    }
                }
            } else {
                quote! {}  // There is no squared norm, so return no impls for norm
            };

            */

            // Implement .normalized() which returns a scaled copy of the object
            // where the bulk norm has been made to equal to 1
            // Also implement .unitized() which returns a scaled copy of the object
            // where the weight norm has been made to 𝟙
            let hat_code = if !obj.is_scalar {
                let type_name = obj.type_name();
                quote! {
                    impl<T: Ring + Recip<Output=T>> Normalized for #type_name
                    where
                        Self: BulkNorm<Output=T>
                    {
                        type Output = Self;
                        fn normalized(self) -> Self {
                            self * self.bulk_norm().recip()
                        }
                    }

                    impl<T: Ring + Recip<Output=T>> Unitized for #type_name
                    where
                        Self: WeightNorm<Output=AntiScalar<T>>
                    {
                        type Output = Self;
                        fn unitized(self) -> Self {
                            self.anti_mul(self.weight_norm().anti_recip())
                        }
                    }
                }
            } else {
                quote! {}
            };

            // Overload +
            let op_trait = quote! { core::ops::Add };
            let op_fn = Ident::new("add", Span::call_site());
            let add_code = gen_binary_operator(
                &basis,
                &objects,
                op_trait,
                op_fn,
                &obj,
                |a, b| generate_symbolic_sum(&basis, a, b, 1, 1),
                true,  // implicit_promotion_to_compound
                None, // alias
            );

            // Overload -
            let op_trait = quote! { core::ops::Sub };
            let op_fn = Ident::new("sub", Span::call_site());
            let sub_code = gen_binary_operator(
                &basis,
                &objects,
                op_trait,
                op_fn,
                &obj,
                |a, b| generate_symbolic_sum(&basis, a, b, 1, -1),
                true,  // implicit_promotion_to_compound
                None, // alias
            );

            // Add a method A.wedge(B) which computes A ∧ B
            let op_trait = quote! { Wedge };
            let op_fn = Ident::new("wedge", Span::call_site());
            let alias_trait = quote! { Join };
            let alias_fn = Ident::new("join", Span::call_site());
            let wedge_product_f = |i: usize, j: usize| {
                let (coef, ix) = geometric_product_multiplication_table[i][j];
                // Select grade s + t
                let s = basis[i].len();
                let t = basis[j].len();
                let u = basis[ix].len();
                let coef = if s + t == u { coef } else { 0 };
                (coef, ix)
            };
            let wedge_product_code = gen_binary_operator(
                &basis,
                &objects,
                op_trait,
                op_fn,
                &obj,
                |a, b| generate_symbolic_product(&basis, a, b, wedge_product_f),
                false, // implicit_promotion_to_compound
                Some((alias_trait, alias_fn)), // alias
            );

            // Add a method A.anti_wedge(B) which computes A ∨ B
            let op_trait = quote! { AntiWedge };
            let op_fn = Ident::new("anti_wedge", Span::call_site());
            let alias_trait = quote! { Meet };
            let alias_fn = Ident::new("meet", Span::call_site());
            let anti_wedge_product_f = |i: usize, j: usize| {

                let (i_coef, i) = right_complement(&right_complement_signs, i);
                let (j_coef, j) = right_complement(&right_complement_signs, j);

                let (prod_coef, ix) = geometric_product_multiplication_table[i][j];
                // Select grade s + t
                let s = basis[i].len();
                let t = basis[j].len();
                let u = basis[ix].len();
                let prod_coef = if s + t == u { prod_coef } else { 0 };

                let (complement_coef, ix) = left_complement(&right_complement_signs, ix);
                (i_coef * j_coef * prod_coef * complement_coef, ix)
            };
            let anti_wedge_product_code = gen_binary_operator(
                &basis,
                &objects,
                op_trait,
                op_fn,
                &obj,
                |a, b| generate_symbolic_product(&basis, a, b, anti_wedge_product_f),
                false, // implicit_promotion_to_compound
                Some((alias_trait, alias_fn)), // alias
            );

            // Add a method A.dot(B) which computes A • B
            let op_trait = quote! { Dot };
            let op_fn = Ident::new("dot", Span::call_site());
            let dot_product = |i: usize, j: usize| (dot_product_multiplication_table[i][j], 0);

            let dot_product_code = gen_binary_operator(
                &basis,
                &objects,
                op_trait,
                op_fn,
                &obj,
                |a, b| generate_symbolic_product(&basis, a, b, dot_product),
                false, // implicit_promotion_to_compound
                None, // alias
            );

            // Add a method A.anti_dot(B) which computes A ∘ B
            let op_trait = quote! { AntiDot };
            let op_fn = Ident::new("anti_dot", Span::call_site());
            let anti_dot_product = |i: usize, j: usize| {
                let (i_coef, i) = right_complement(&right_complement_signs, i);
                let (j_coef, j) = right_complement(&right_complement_signs, j);

                let (prod_coef, ix) = dot_product(i, j);

                let (complement_coef, ix) = left_complement(&right_complement_signs, ix);
                (i_coef * j_coef * prod_coef * complement_coef, ix)
            };

            let anti_dot_product_code = gen_binary_operator(
                &basis,
                &objects,
                op_trait,
                op_fn,
                &obj,
                |a, b| generate_symbolic_product(&basis, a, b, anti_dot_product),
                false, // implicit_promotion_to_compound
                None, // alias
            );

            // Implement the geometric product ⟑
            let op_trait = quote! { WedgeDot };
            let op_fn = Ident::new("wedge_dot", Span::call_site());

            let wedge_dot_product_code = gen_binary_operator(
                &basis,
                &objects,
                op_trait,
                op_fn,
                &obj,
                |a, b| generate_symbolic_product(&basis, a, b, |i: usize, j: usize| geometric_product_multiplication_table[i][j]),
                false, // implicit_promotion_to_compound
                None, // alias
            );

            // Implement the geometric antiproduct ⟇
            let op_trait = quote! { AntiWedgeDot };
            let op_fn = Ident::new("anti_wedge_dot", Span::call_site());
            let alias_trait = quote! { Compose };
            let alias_fn = Ident::new("compose", Span::call_site());
            let anti_wedge_dot_product = |i: usize, j: usize| geometric_anti_product_multiplication_table[i][j];

            let anti_wedge_dot_product_code = gen_binary_operator(
                &basis,
                &objects,
                op_trait,
                op_fn,
                &obj,
                |a, b| generate_symbolic_product(&basis, a, b, anti_wedge_dot_product),
                false, // implicit_promotion_to_compound
                Some((alias_trait, alias_fn)), // alias
            );

            // Overload * for scalar multiplication
            let op_trait = quote! { core::ops::Mul };
            let op_fn = Ident::new("mul", Span::call_site());
            let scalar_product = |i: usize, j: usize| {
                let (coef, ix) = geometric_product_multiplication_table[i][j];
                // Require s or t be a scalar
                let s = basis[i].len();
                let t = basis[j].len();
                let coef = if s == 0 || t == 0 { coef } else { 0 };
                (coef, ix)
            };
            let scalar_product_code = gen_binary_operator(
                &basis,
                &objects,
                op_trait,
                op_fn,
                &obj,
                |a, b| generate_symbolic_product(&basis, a, b, scalar_product),
                true, // implicit_promotion_to_compound
                None, // alias
            );

            // Implement A.anti_mul(B) for anti-scalar multiplication
            let op_trait = quote! { AntiMul };
            let op_fn = Ident::new("anti_mul", Span::call_site());
            let anti_scalar_f = |i: usize, j: usize| {
                let (i_coef, i) = right_complement(&right_complement_signs, i);
                let (j_coef, j) = right_complement(&right_complement_signs, j);

                let (prod_coef, ix) = scalar_product(i, j);

                let (complement_coef, ix) = left_complement(&right_complement_signs, ix);
                (i_coef * j_coef * prod_coef * complement_coef, ix)
            };
            let anti_scalar_product_code = gen_binary_operator(
                &basis,
                &objects,
                op_trait,
                op_fn,
                &obj,
                |a, b| generate_symbolic_product(&basis, a, b, anti_scalar_f),
                true, // implicit_promotion_to_compound
                None, // alias
            );

            // Add a method A.bulk_expansion(B) which computes A ∧ B★
            let op_trait = quote! { BulkExpansion };
            let op_fn = Ident::new("bulk_expansion", Span::call_site());
            let bulk_expansion_f = |i: usize, j: usize| {
                let (coef_dual, j) = bulk_dual_f(j);
                let (coef_prod, ix) = wedge_product_f(i, j);
                (coef_dual * coef_prod, ix)
            };
            let bulk_expansion_code = gen_binary_operator(
                &basis,
                &objects,
                op_trait,
                op_fn,
                &obj,
                |a, b| generate_symbolic_product(&basis, a, b, bulk_expansion_f),
                false, // implicit_promotion_to_compound
                None, // alias
            );

            // Add a method A.weight_expansion(B) which computes A ∧ B☆
            let op_trait = quote! { WeightExpansion };
            let op_fn = Ident::new("weight_expansion", Span::call_site());
            let alias_trait = quote! { SupersetOrthogonalTo };
            let alias_fn = Ident::new("superset_orthogonal_to", Span::call_site());
            let weight_expansion_f = |i: usize, j: usize| {
                let (coef_dual, j) = weight_dual_f(j);
                let (coef_prod, ix) = wedge_product_f(i, j);
                (coef_dual * coef_prod, ix)
            };
            let weight_expansion_code = gen_binary_operator(
                &basis,
                &objects,
                op_trait,
                op_fn,
                &obj,
                |a, b| generate_symbolic_product(&basis, a, b, weight_expansion_f),
                false, // implicit_promotion_to_compound
                Some((alias_trait, alias_fn)), // alias
            );

            // Add a method A.bulk_contraction(B) which computes A ∨ B★
            let op_trait = quote! { BulkContraction };
            let op_fn = Ident::new("bulk_contraction", Span::call_site());
            let bulk_contraction_f = |i: usize, j: usize| {
                let (coef_dual, j) = bulk_dual_f(j);
                let (coef_prod, ix) = anti_wedge_product_f(i, j);
                (coef_dual * coef_prod, ix)
            };
            let bulk_contraction_code = gen_binary_operator(
                &basis,
                &objects,
                op_trait,
                op_fn,
                &obj,
                |a, b| generate_symbolic_product(&basis, a, b, bulk_contraction_f),
                false, // implicit_promotion_to_compound
                None, // alias
            );

            // Add a method A.weight_contraction(B) which computes A ∨ B☆
            let op_trait = quote! { WeightContraction };
            let op_fn = Ident::new("weight_contraction", Span::call_site());
            let alias_trait = quote! { SubsetOrthogonalTo };
            let alias_fn = Ident::new("subset_orthogonal_to", Span::call_site());
            let weight_contraction_f = |i: usize, j: usize| {
                let (coef_dual, j) = weight_dual_f(j);
                let (coef_prod, ix) = anti_wedge_product_f(i, j);
                (coef_dual * coef_prod, ix)
            };
            let weight_contraction_code = gen_binary_operator(
                &basis,
                &objects,
                op_trait,
                op_fn,
                &obj,
                |a, b| generate_symbolic_product(&basis, a, b, weight_contraction_f),
                false, // implicit_promotion_to_compound
                Some((alias_trait, alias_fn)), // alias
            );

            // Add a method A.anti_commutator(B) which computes (A ⟇ B - B ⟇ A) / 2
            let op_trait = quote! { AntiCommutator };
            let op_fn = Ident::new("anti_commutator", Span::call_site());
            let commutator_product = |i: usize, j: usize| {
                let (coef_1, ix) = geometric_anti_product_multiplication_table[i][j];
                let (coef_2, ix_2) = geometric_anti_product_multiplication_table[j][i];

                let coef = coef_1 - coef_2;

                assert!(ix == ix_2);
                assert!(coef % 2 == 0);
                let coef = coef / 2;

                (coef, ix)
            };
            let commutator_product_code = gen_binary_operator(
                &basis,
                &objects,
                op_trait,
                op_fn,
                &obj,
                |a, b| generate_symbolic_product(&basis, a, b, commutator_product),
                false, // implicit_promotion_to_compound
                None, // alias
            );

            // Add a method A.transform(B) which computes B̰ ⟇ A ⟇ B
            let op_trait = quote! { Transform };
            let op_fn = Ident::new("transform", Span::call_site());
            let transform_1 = |i: usize, j: usize| {
                // Compute first half of B̰ ⟇ A ⟇ B
                // Where i maps to B, and j maps to A.
                // In part 2, we will compute the geometric antiproduct of this intermediate result
                // with B

                let (coef_rev, i) = anti_reverse_f(i);
                let (coef_prod, ix_result) = geometric_anti_product_multiplication_table[i][j];
                (coef_rev * coef_prod, ix_result)
            };
            let transform_2 = |i: usize, j: usize| {
                // Compute second half of B̰ ⟇ A ⟇ B
                // In part 1, we computed the intermediate result B̰ ⟇ A which maps to i here.
                // j maps to B.

                geometric_anti_product_multiplication_table[i][j]
            };
            let transform_code = gen_binary_operator(
                &basis,
                &objects,
                op_trait,
                op_fn,
                &obj,
                |a, b| {
                    generate_symbolic_double_product(
                        &basis,
                        a,
                        b,
                        transform_1,
                        transform_2,
                    )
                },
                false, // implicit_promotion_to_compound
                None, // alias
            );

            // Implement A.projection(B) which computes B ∨ (A ∧  B☆)
            let op_trait = quote! { Projection };
            let op_fn = Ident::new("projection", Span::call_site());
            let projection_product_1 = |i: usize, j: usize| {
                // Compute second half of B ∨ (A ∧ B☆)
                // Where i maps to B, and j maps to A.
                // In part 2, we will compute the geometric product of B with this intermediate result
                weight_expansion_f(j, i)
            };
            let projection_product_2 = |i: usize, j: usize| {
                // Compute second half of B ∨ (A ∧ B☆)
                // In part 1, we computed the intermediate result A ∧ B☆ which maps to i here.
                // j maps to B.
                anti_wedge_product_f(j, i)
            };
            let projection_code = gen_binary_operator(
                &basis,
                &objects,
                op_trait,
                op_fn,
                &obj,
                |a, b| {
                    generate_symbolic_double_product(
                        &basis,
                        a,
                        b,
                        projection_product_1,
                        projection_product_2,
                    )
                },
                false, // implicit_promotion_to_compound
                None, // alias
            );

            // Implement A.anti_projection(B) which computes B ∧ (A ∨ B☆)
            let op_trait = quote! { AntiProjection };
            let op_fn = Ident::new("anti_projection", Span::call_site());
            let anti_projection_product_1 = |i: usize, j: usize| {
                // Compute second half of B ∧ (A ∨ B☆)
                // Where i maps to B, and j maps to A.
                // In part 2, we will compute the geometric product of B with this intermediate result
                weight_contraction_f(j, i)
            };
            let anti_projection_product_2 = |i: usize, j: usize| {
                // Compute second half of B ∧ (A ∨ B☆)
                // In part 1, we computed the intermediate result A ∨ B☆ which maps to i here.
                // j maps to B.
                wedge_product_f(j, i)
            };
            let anti_projection_code = gen_binary_operator(
                &basis,
                &objects,
                op_trait,
                op_fn,
                &obj,
                |a, b| {
                    generate_symbolic_double_product(
                        &basis,
                        a,
                        b,
                        anti_projection_product_1,
                        anti_projection_product_2,
                    )
                },
                false, // implicit_promotion_to_compound
                None, // alias
            );

            // Implement A.central_projection(B) which computes B ∨ (A ∧ B★)
            let op_trait = quote! { CentralProjection };
            let op_fn = Ident::new("central_projection", Span::call_site());
            let central_projection_product_1 = |i: usize, j: usize| {
                // Compute second half of B ∨ (A ∧ B★)
                // Where i maps to B, and j maps to A.
                // In part 2, we will compute the geometric product of B with this intermediate result
                bulk_expansion_f(j, i)
            };
            let central_projection_product_2 = |i: usize, j: usize| {
                // Compute second half of B ∨ (A ∧ B★)
                // In part 1, we computed the intermediate result A ∧ B★ which maps to i here.
                // j maps to B.
                anti_wedge_product_f(j, i)
            };
            let central_projection_code = gen_binary_operator(
                &basis,
                &objects,
                op_trait,
                op_fn,
                &obj,
                |a, b| {
                    generate_symbolic_double_product(
                        &basis,
                        a,
                        b,
                        central_projection_product_1,
                        central_projection_product_2,
                    )
                },
                false, // implicit_promotion_to_compound
                None, // alias
            );

            // Implement A.central_anti_projection(B) which computes B ∧ (A ∨ B★)
            let op_trait = quote! { CentralAntiProjection };
            let op_fn = Ident::new("central_anti_projection", Span::call_site());
            let central_anti_projection_product_1 = |i: usize, j: usize| {
                // Compute second half of B ∧ (A ∨ B★)
                // Where i maps to B, and j maps to A.
                // In part 2, we will compute the geometric product of B with this intermediate result
                bulk_contraction_f(j, i)
            };
            let central_anti_projection_product_2 = |i: usize, j: usize| {
                // Compute second half of B ∧ (A ∨ B★)
                // In part 1, we computed the intermediate result A ∨ B★ which maps to i here.
                // j maps to B.
                wedge_product_f(j, i)
            };
            let central_anti_projection_code = gen_binary_operator(
                &basis,
                &objects,
                op_trait,
                op_fn,
                &obj,
                |a, b| {
                    generate_symbolic_double_product(
                        &basis,
                        a,
                        b,
                        central_anti_projection_product_1,
                        central_anti_projection_product_2,
                    )
                },
                false, // implicit_promotion_to_compound
                None, // alias
            );

            // Implement motor_to which computes A̰ ⟇ B
            let op_trait = quote! { MotorTo };
            let op_fn = Ident::new("motor_to", Span::call_site());
            let motor_to_f = |i: usize, j: usize| {
                let (coef_rev, i) = anti_reverse_f(i);
                let (coef_prod, ix_result) = geometric_anti_product_multiplication_table[i][j];
                (coef_rev * coef_prod, ix_result)
            };

            let motor_to_code = gen_binary_operator(
                &basis,
                &objects,
                op_trait,
                op_fn,
                &obj,
                |a, b| generate_symbolic_product(&basis, a, b, motor_to_f),
                true, // implicit_promotion_to_compound
                None, // alias
            );

            quote! {
                // ===========================================================================
                // #name
                // ===========================================================================

                #clone_code
                #copy_code
                #debug_code
                #display_code
                #partial_eq_code
                #eq_code
                #neg_code
                #reverse_code
                #anti_reverse_code
                #bulk_code
                #weight_code
                #bulk_dual_code
                #weight_dual_code
                #right_complement_code
                #left_complement_code
                #bulk_norm_squared_code
                #bulk_norm_code
                #weight_norm_squared_code
                #weight_norm_code
                #hat_code
                #add_code
                #sub_code
                #wedge_product_code
                #anti_wedge_product_code
                #dot_product_code
                #anti_dot_product_code
                #wedge_dot_product_code
                #anti_wedge_dot_product_code
                #scalar_product_code
                #anti_scalar_product_code
                #commutator_product_code
                #bulk_expansion_code
                #weight_expansion_code
                #bulk_contraction_code
                #weight_contraction_code
                #projection_code
                #anti_projection_code
                #central_projection_code
                #central_anti_projection_code
                #transform_code
                #motor_to_code
            }
        })
        .collect();

    quote! {
        #struct_code
        #anti_scalar_ops
        #impl_code
    }
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
