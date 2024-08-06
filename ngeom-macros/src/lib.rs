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
                            tokens.append_all(quote! { T::two() * });
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

    let basis_element_count = basis.len();

    fn coefficient_ident(var: &str, basis: &[usize]) -> Ident {
        let basis_pretty: String = basis.iter().map(|e| format!("{}", e)).collect();
        Ident::new(&format!("{}{}", var, basis_pretty), Span::call_site())
    }

    // We will represent a multivector as an array of coefficients on the basis elements.
    // e.g. in 2D, there are 1 + 3 + 3 + 1 = 8 basis elements,
    // and a full multivector uses all of them: [1, 1, 1, 1, 1, 1, 1, 1]
    // An object such as a Bivector would only use a few of them: [0, 0, 0, 0, 1, 1, 1, 0]

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

    // Add even & odd sub-algebras to the object set
    for parity in 0..2 {
        let select_components: Vec<_> = basis.iter().map(|b| b.len() % 2 == parity).collect();
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

    fn gen_unary_operator<F: Fn(usize) -> (isize, usize)>(
        basis: &Vec<Vec<usize>>,
        objects: &[Object],
        op_trait: TokenStream,
        op_fn: Ident,
        obj: &Object,
        op: F,
        output_self: bool,
    ) -> TokenStream {
        if obj.is_scalar {
            // Do not generate operations with the scalar being the LHS--
            // typically because these would violate rust's orphan rule
            // or result in conflicting trait implementations
            return quote! {};
        }

        // Generate the vec of symbolic expressions (one for each basis element coefficient)
        // that computes the given unary operator
        let mut expressions: Vec<SymbolicSumExpr> =
            vec![Default::default(); obj.select_components.len()];
        for i in obj
            .select_components
            .iter()
            .enumerate()
            .filter_map(|(i, is_selected)| is_selected.then_some(i))
        {
            let (coef, result_basis_ix) = op(i);

            expressions[result_basis_ix] +=
                SymbolicProdExpr(coef, vec![Symbol(coefficient_ident("a", &basis[i]))]);
        }

        let expressions: Vec<_> = expressions
            .into_iter()
            .map(|expr| expr.simplify())
            .collect();

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

        let output_type = if output_self {
            quote! {}
        } else {
            quote! { type Output = #output_type_name; }
        };

        quote! {
            impl < T: Ring > #op_trait for #type_name {
                #output_type

                fn #op_fn (self) -> #output_type_name {
                    #return_expr
                }
            }
        }
    }

    // Generate multiplication table for the basis elements
    // Each entry is a tuple of the (coefficient, basis_index)
    // e.g. (1, 0) means the multiplication result is 1 * scalar = 1
    let multiplication_table: Vec<Vec<(isize, usize)>> = {
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
    // reflect, etc.
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
        output_self: bool,
        implicit_promotion_to_compound: bool,
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

                if output_object.is_scalar && expressions[0].0.len() == 0 {
                    // This operation unconditionally returns 0,
                    // so invoking it is probably a type error--
                    // do not generate code for it
                    return quote! {};
                }

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

                let output_type_name = &output_object.type_name_colons();
                let rhs_type_name = &rhs_obj.type_name();
                let lhs_type_name = &lhs_obj.type_name();

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

                let output_type = if output_self {
                    quote! {}
                } else {
                    quote! { type Output = #output_type_name; }
                };

                quote! {
                    impl < T: Ring > #op_trait < #rhs_type_name >  for #lhs_type_name {
                        #output_type

                        fn #op_fn (self, r: #rhs_type_name) -> #output_type_name {
                            #return_expr
                        }
                    }
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
                        #[derive(Clone, Copy)]
                        pub struct #name <T: Ring> {
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

    let impl_code: TokenStream = objects
        .iter()
        .map(|obj| {
            // Overload unary -
            let op_trait = quote! { std::ops::Neg };
            let op_fn = Ident::new("neg", Span::call_site());
            let neg_f = |i: usize| (-1, i);
            let neg_code =
                gen_unary_operator(&basis, &objects, op_trait, op_fn, &obj, neg_f, false);

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
            let reverse_code =
                gen_unary_operator(&basis, &objects, op_trait, op_fn, &obj, reverse_f, true);

            // Add a method A.dual()
            let op_trait = quote! { Dual };
            let op_fn = Ident::new("dual", Span::call_site());
            let dual_f = |i: usize| {
                // We have set up our basis carefully to allow this
                (1, basis.len() - i - 1)
            };
            let dual_code =
                gen_unary_operator(&basis, &objects, op_trait, op_fn, &obj, dual_f, false);

            // Add methods: norm_squared() and inorm_squared()
            // as well as norm() and inorm() if the scalar type implements sqrt()
            let norm_code = {
                if obj.is_scalar {
                    // Do not generate norm operations for the scalar type
                    // since this will result in conflicting trait implementations
                    return quote! {};
                }

                let gen_norm_squared_return_expr = |dual| {
                    let mut expressions: Vec<SymbolicSumExpr> =
                        vec![Default::default(); obj.select_components.len()];
                    // The squared norms are computed like a product
                    // so this code looks similar to the generic product code
                    for i in obj
                        .select_components
                        .iter()
                        .enumerate()
                        .filter_map(|(i, is_selected)| is_selected.then_some(i))
                    {
                        for j in obj
                            .select_components
                            .iter()
                            .enumerate()
                            .filter_map(|(j, is_selected)| is_selected.then_some(j))
                        {
                            // Compute the product A * ~A
                            // or dual(A) * ~dual(A) depending on dual flag

                            let symbols = vec![
                                Symbol(coefficient_ident("a", &basis[i])),
                                Symbol(coefficient_ident("a", &basis[j])),
                            ];

                            let (i, j) = if dual {
                                (basis.len() - i - 1, basis.len() - j - 1)
                            } else {
                                (i, j)
                            };

                            let reverse_coef = match (basis[j].len() / 2) % 2 {
                                0 => 1,
                                1 => -1,
                                _ => panic!("Expected parity to be 0 or 1"),
                            };

                            let (coef, ix_result) = multiplication_table[i][j];

                            expressions[ix_result] +=
                                SymbolicProdExpr(coef * reverse_coef, symbols);
                        }
                    }

                    let expressions: Vec<_> = expressions
                        .into_iter()
                        .map(|expr| expr.simplify())
                        .collect();

                    // Ensure the result is a non-zero scalar
                    // TODO: should it be?
                    //assert!(
                    //    expressions
                    //        .iter()
                    //        .enumerate()
                    //        .all(|(i, SymbolicSumExpr(terms))| if i == 0 {
                    //            true
                    //        } else {
                    //            if terms.len() != 0 {
                    //                let expr = SymbolicSumExpr(terms.clone());
                    //                println!("Obj {} norm {}: Coef {:?} was {}", obj.name, dual, i, quote! { #expr });
                    //            }
                    //            terms.len() == 0
                    //        }),
                    //    "Expression for squared norm was not a non-zero scalar"
                    //);

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

                    // Take expressions[0]
                    expressions.into_iter().next().unwrap()
                };

                let type_name = obj.type_name_colons();
                let norm_squared_return_expr = gen_norm_squared_return_expr(false);
                let inorm_squared_return_expr = gen_norm_squared_return_expr(true);

                quote! {
                    impl < T: Ring > NormSquared for #type_name {
                        type Output = T;

                        fn norm_squared (self) -> T {
                            #norm_squared_return_expr
                        }
                    }

                    impl < T: Ring > INormSquared for #type_name {
                        type Output = T;
                        fn inorm_squared (self) -> T {
                            #inorm_squared_return_expr
                        }
                    }
                    // TODO: Some norms can be optimized to avoid needing Sqrt
                    impl < T: Ring + Sqrt > Norm for #type_name {}
                    impl < T: Ring + Sqrt > INorm for #type_name {}
                }
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
                false, // output_self
                true,  // implicit_promotion_to_compound
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
                false, // output_self
                true,  // implicit_promotion_to_compound
            );

            // Overload ^ with the wedge product
            let op_trait = quote! { core::ops::BitXor };
            let op_fn = Ident::new("bitxor", Span::call_site());
            let wedge_product = |i: usize, j: usize| {
                let (coef, ix) = multiplication_table[i][j];
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
                |a, b| generate_symbolic_product(&basis, a, b, wedge_product),
                false, // output_self
                false, // implicit_promotion_to_compound
            );

            // Overload & with the vee product
            let op_trait = quote! { core::ops::BitAnd };
            let op_fn = Ident::new("bitand", Span::call_site());
            let vee_product = |i: usize, j: usize| {
                let dual = |ix| basis.len() - ix - 1;

                let (coef, ix) = multiplication_table[dual(i)][dual(j)];
                // Select grade s + t
                let s = basis[i].len();
                let t = basis[j].len();
                let u = basis[ix].len();
                let coef = if s + t == u { coef } else { 0 };

                let ix = dual(ix);
                (coef, ix)
            };
            let vee_product_code = gen_binary_operator(
                &basis,
                &objects,
                op_trait,
                op_fn,
                &obj,
                |a, b| generate_symbolic_product(&basis, a, b, vee_product),
                false, // output_self
                false, // implicit_promotion_to_compound
            );

            // Overload | with the dot product
            let op_trait = quote! { core::ops::BitOr };
            let op_fn = Ident::new("bitor", Span::call_site());
            let dot_product = |i: usize, j: usize| {
                let (coef, ix) = multiplication_table[i][j];
                // Select grade |s - t|
                let s = basis[i].len();
                let t = basis[j].len();
                let u = basis[ix].len();
                let coef = if s + u == t || t + u == s { coef } else { 0 };
                (coef, ix)
            };

            let dot_product_code = gen_binary_operator(
                &basis,
                &objects,
                op_trait,
                op_fn,
                &obj,
                |a, b| generate_symbolic_product(&basis, a, b, dot_product),
                false, // output_self
                false, // implicit_promotion_to_compound
            );

            // Overload * with the geometric product
            let op_trait = quote! { core::ops::Mul };
            let op_fn = Ident::new("mul", Span::call_site());
            let geometric_product = |i: usize, j: usize| multiplication_table[i][j];
            let geometric_product_code = gen_binary_operator(
                &basis,
                &objects,
                op_trait,
                op_fn,
                &obj,
                |a, b| generate_symbolic_product(&basis, a, b, geometric_product),
                false, // output_self
                false, // implicit_promotion_to_compound
            );

            // Add a method A.cross(B) which computes (A * B - B * A) / 2
            let op_trait = quote! { Commutator };
            let op_fn = Ident::new("cross", Span::call_site());
            let commutator_product = |i: usize, j: usize| {
                let (coef_1, ix) = multiplication_table[i][j];
                let (coef_2, ix_2) = multiplication_table[j][i];

                let coef = coef_1 + coef_2;

                assert!(ix == ix_2);
                assert!(coef % 2 == 0);

                (coef / 2, ix)
            };
            let commutator_product_code = gen_binary_operator(
                &basis,
                &objects,
                op_trait,
                op_fn,
                &obj,
                |a, b| generate_symbolic_product(&basis, a, b, commutator_product),
                false, // output_self
                false, // implicit_promotion_to_compound
            );

            // Add a method A.transform(B) which computes B * A * ~B
            let op_trait = quote! { Transform };
            let op_fn = Ident::new("transform", Span::call_site());
            let transform_product_1 = |i: usize, j: usize| {
                // Compute first half of B * A * ~B
                // Where i maps to B, and j maps to A.
                // In part 2, we will compute the geometric product of this intermediate result
                // with ~B

                multiplication_table[i][j]
            };
            let transform_product_2 = |i: usize, j: usize| {
                // Compute second half of B * A * ~B
                // In part 1, we computed the intermediate result B * A which maps to i here.
                // j maps to B.
                let reverse_coef = |ix: usize| match (basis[ix].len() / 2) % 2 {
                    0 => 1,
                    1 => -1,
                    _ => panic!("Expected parity to be 0 or 1"),
                };

                let (coef, ix_result) = multiplication_table[i][j];

                (coef * reverse_coef(j), ix_result)
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
                        transform_product_1,
                        transform_product_2,
                    )
                },
                true,  // output_self
                false, // implicit_promotion_to_compound
            );

            // Add a method A.reflect(B) which computes B * A * B
            let op_trait = quote! { Reflect };
            let op_fn = Ident::new("reflect", Span::call_site());
            let reflect_product_1 = |i: usize, j: usize| {
                // Compute first half of B * A * B
                // Where i maps to B, and j maps to A.
                // In part 2, we will compute the geometric product of this intermediate result
                // with B

                multiplication_table[i][j]
            };
            let reflect_product_2 = |i: usize, j: usize| {
                // Compute second half of B * A * B
                // In part 1, we computed the intermediate result B * A which maps to i here.
                // j maps to B.
                multiplication_table[i][j]
            };
            let reflect_code = gen_binary_operator(
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
                        reflect_product_1,
                        reflect_product_2,
                    )
                },
                true,  // output_self
                false, // implicit_promotion_to_compound
            );

            // Add a method A.project(B) which computes (B . A) * B
            let op_trait = quote! { Project };
            let op_fn = Ident::new("project", Span::call_site());
            let project_product_1 = |i: usize, j: usize| {
                // Compute first half of (B . A) * B
                // Where i maps to B, and j maps to A.
                // In part 2, we will compute the geometric product of this intermediate result
                // with B

                let (coef, ix) = multiplication_table[i][j];
                // Select grade |s - t|
                let s = basis[i].len();
                let t = basis[j].len();
                let u = basis[ix].len();
                let coef = if s + u == t || t + u == s { coef } else { 0 };
                (coef, ix)
            };
            let project_product_2 = |i: usize, j: usize| {
                // Compute second half of (B . A) * B
                // In part 1, we computed the intermediate result B . A which maps to i here.
                // j maps to B.
                multiplication_table[i][j]
            };
            let project_code = gen_binary_operator(
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
                        project_product_1,
                        project_product_2,
                    )
                },
                true,  // output_self
                false, // implicit_promotion_to_compound
            );

            quote! {
                // ===========================================================================
                // #name
                // ===========================================================================

                // TODO: NormSquared & Norm
                // TODO: exp
                #neg_code
                #reverse_code
                #dual_code
                #norm_code
                #add_code
                #sub_code
                #wedge_product_code
                #vee_product_code
                #dot_product_code
                #geometric_product_code
                #commutator_product_code
                #project_code
                #reflect_code
                #transform_code
            }
        })
        .collect();

    quote! {
        #struct_code
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
