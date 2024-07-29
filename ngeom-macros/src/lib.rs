extern crate proc_macro;
use core::cmp::Ordering;
use core::iter::Sum;
use core::ops::{Add, Mul};
use nalgebra::DVector;
use proc_macro2::Span;
use proc_macro2::TokenStream;
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
            self.0.to_tokens(tokens)
        }
    }

    impl PartialOrd for SymbolicProdExpr {
        fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
            Some(self.cmp(other))
        }
    }

    impl Ord for SymbolicProdExpr {
        fn cmp(&self, other: &Self) -> Ordering {
            self.1.cmp(&other.1).then_with(|| self.0.cmp(&other.0))
        }
    }

    impl Mul<SymbolicProdExpr> for SymbolicProdExpr {
        type Output = SymbolicProdExpr;
        fn mul(mut self, mut r: SymbolicProdExpr) -> SymbolicProdExpr {
            self.0 *= r.0;
            self.1.append(&mut r.1);
            self
        }
    }

    impl SymbolicProdExpr {
        fn simplify(mut self) -> Self {
            // Sort expression
            if self.0 == 0 {
                self.1.clear();
            } else {
                self.1.sort();
            }
            self
        }
    }

    impl ToTokens for SymbolicSumExpr {
        fn to_tokens(&self, tokens: &mut TokenStream) {
            if self.0.len() == 0 {
                tokens.append_all(quote! { ScalarRing::zero() });
            } else {
                for (count, prod_expr) in self.0.iter().enumerate() {
                    let SymbolicProdExpr(coef, prod_terms) = prod_expr;
                    let coef = *coef;

                    if coef >= 0 {
                        if count != 0 {
                            tokens.append_all(quote! { + });
                        }
                    } else {
                        tokens.append_all(quote! { - });
                    }

                    if prod_terms.len() == 0 {
                        // If there are no symbols in the product, then this is a scalar
                        if coef == 0 {
                            tokens.append_all(quote! { ScalarRing::zero() });
                        } else if coef.abs() == 1 {
                            tokens.append_all(quote! {
                                ScalarRing::one()
                            });
                        } else {
                            panic!("Scalar was not 0, -1 or 1");
                        }
                    } else {
                        // There are symbols in the product
                        if coef.abs() != 1 {
                            panic!("Coefficient on symbols was not -1 or 1");
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

    //impl From<isize> for SymbolicProdExpr {
    //    fn from(x: isize) -> SymbolicProdExpr {
    //        SymbolicProdExpr(x, vec![])
    //    }
    //}

    impl SymbolicSumExpr {
        fn simplify(mut self) -> Self {
            // Simplify all products
            self.0 = self.0.into_iter().map(|prod| prod.simplify()).collect();

            // Sort expression by symbolic values
            self.0.sort();

            // Combine adjacent terms whose symbolic parts are equal
            let mut new_expression = vec![];
            let mut prev_coef = 0;
            let mut prev_symbols = vec![];
            for SymbolicProdExpr(coef, symbols) in self.0.into_iter() {
                if prev_symbols == symbols {
                    prev_coef += coef;
                } else {
                    new_expression.push(SymbolicProdExpr(prev_coef, prev_symbols));
                    prev_coef = coef;
                    prev_symbols = symbols;
                }
            }
            new_expression.push(SymbolicProdExpr(prev_coef, prev_symbols));

            self.0 = new_expression;

            // Remove all products with coefficient = 0
            self.0.retain(|prod| prod.0 != 0);

            self
        }
    }

    //impl Add<SymbolicSumExpr> for SymbolicSumExpr {
    //    type Output = SymbolicSumExpr;
    //    fn add(mut self, mut r: SymbolicSumExpr) -> SymbolicSumExpr {
    //        self.0.append(&mut r.0);
    //        self
    //    }
    //}

    impl Mul<SymbolicSumExpr> for SymbolicSumExpr {
        type Output = SymbolicSumExpr;
        fn mul(self, r: SymbolicSumExpr) -> SymbolicSumExpr {
            let SymbolicSumExpr(l) = self;
            SymbolicSumExpr(
                l.iter()
                    .flat_map(|lp| r.0.iter().map(|rp| lp.clone() * rp.clone()))
                    .collect(),
            )
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

    fn coefficient_symbol(var: &str, basis: &[usize]) -> Symbol {
        let basis_pretty: String = basis.iter().map(|e| format!("{}", e)).collect();
        Symbol(Ident::new(
            &format!("{}{}", var, basis_pretty),
            Span::call_site(),
        ))
    }

    // We will represent a multivector as an array of coefficients on the basis elements.
    // e.g. in 2D, there are 1 + 3 + 3 + 1 = 8 basis elements,
    // and a full multivector uses all of them: [1, 1, 1, 1, 1, 1, 1, 1]
    // An object such as a Bivector would only use a few of them: [0, 0, 0, 0, 1, 1, 1, 0]

    struct Object {
        name: Ident,
        select_components: DVector<bool>,
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
        let select_components =
            DVector::<bool>::from_iterator(basis_element_count, basis.iter().map(|b| b.len() == k));
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
        let select_components = DVector::<bool>::from_iterator(
            basis_element_count,
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
            is_scalar: false,
            is_compound: true,
        });
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

    // Function to generate a product of two objects, e.g. geometric product, wedge product, etc.
    // The result will be a list of symbolic expressions,
    // one for each coefficient value in the resultant multivector.
    fn generate_symbolic_product<F: Fn(usize, usize) -> (isize, usize)>(
        basis: &Vec<Vec<usize>>,
        select_components_a: &DVector<bool>,
        select_components_b: &DVector<bool>,
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
                if coef != 0 {
                    result[result_basis_ix].0.push(SymbolicProdExpr(
                        coef,
                        vec![
                            coefficient_symbol("a", &basis[i]),
                            coefficient_symbol("b", &basis[j]),
                        ],
                    ));
                }
            }
        }

        result.into_iter().map(|expr| expr.simplify()).collect()
    }

    //fn simplify_symbolic_product(
    //    mut product: Vec<Vec<(isize, usize, usize)>>,
    //) -> Vec<Vec<(isize, usize, usize)>> {
    //    // Simpify each expression by combining repeated entries
    //    for expression in product.iter_mut() {
    //        // Sort expression so duplicates are adjacent
    //        expression.sort_by_key(|&(_, i, j)| (i, j));

    //        // Combine duplicates and remove any with coefficient=0
    //        let mut new_expression = vec![];
    //        let mut prev_ij = None;
    //        let mut prev_coef = 0;
    //        for &(coef, i, j) in expression.iter() {
    //            if Some((i, j)) == prev_ij {
    //                prev_coef += coef;
    //            } else {
    //                if let Some((prev_i, prev_j)) = prev_ij {
    //                    if prev_coef != 0 {
    //                        new_expression.push((prev_coef, prev_i, prev_j));
    //                    }
    //                }
    //                prev_coef = coef;
    //                prev_ij = Some((i, j));
    //            }
    //        }
    //        if let Some((prev_i, prev_j)) = prev_ij {
    //            if prev_coef != 0 {
    //                new_expression.push((prev_coef, prev_i, prev_j));
    //            }
    //        }

    //        *expression = new_expression;
    //    }
    //    product
    //}

    fn gen_binary_operator<F: Fn(&DVector<bool>, &DVector<bool>) -> Vec<SymbolicSumExpr>>(
        basis: &Vec<Vec<usize>>,
        objects: &[Object],
        op_trait: TokenStream,
        op_fn: Ident,
        lhs_obj: &Object,
        op: F,
        output_self: bool,
    ) -> TokenStream {
        objects
            .iter()
            .map(|rhs_obj| {
                if lhs_obj.is_scalar {
                    // Do not generate operations with the scalar being the LHS--
                    // these violate rust's orphan rule
                    return quote! {};
                }

                let expressions = op(&lhs_obj.select_components, &rhs_obj.select_components);

                let select_output_components = DVector::from_iterator(
                    expressions.len(),
                    expressions.iter().map(|e| e.0.len() != 0),
                );
                let output_object = objects.iter().find(|o| {
                    o.select_components
                        .iter()
                        .zip(select_output_components.into_iter())
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

                if output_object.is_compound
                    && !(lhs_obj.is_scalar
                        || rhs_obj.is_scalar
                        || (lhs_obj.is_compound && rhs_obj.is_compound))
                {
                    // Do not create compound objects unintentionally. Only allow:
                    // 1. Products of compound objects and scalars
                    // 2. Products of compound objects and other compound objects
                    return quote! {};
                }

                let output_type_name = &output_object.type_name_colons();
                let rhs_type_name = &rhs_obj.type_name();
                let lhs_type_name = &lhs_obj.type_name();

                //let code_for_expression = |expression: &Vec<(isize, usize, usize)>| -> TokenStream {
                //    if expression.len() == 0 {
                //        return quote! { T::zero() };
                //    }
                //    expression
                //        .iter()
                //        .enumerate()
                //        .map(|(count, &(coef, i, j))| {
                //            let i_expr = if lhs_obj.is_scalar {
                //                assert!(i == 0);
                //                quote! {
                //                    self
                //                }
                //            } else {
                //                let ident = coefficient_identifier("a", &basis[i]);
                //                quote! {
                //                    self.#ident
                //                }
                //            };
                //            let j_expr = if rhs_obj.is_scalar {
                //                assert!(j == 0);
                //                quote! {
                //                    r
                //                }
                //            } else {
                //                let ident = coefficient_identifier("a", &basis[j]);
                //                quote! {
                //                    r.#ident
                //                }
                //            };
                //            let plus = if count == 0 {
                //                quote! {}
                //            } else {
                //                quote! {+}
                //            };
                //            let minus = quote! {-};
                //            if coef == 1 {
                //                quote! { #plus #i_expr * #j_expr}
                //            } else if coef == -1 {
                //                quote! { #minus #i_expr * #j_expr}
                //            } else if coef > 0 {
                //                panic!("Non-unity coefficient");
                //                //let coef_lit = LitInt::new(&format!("{}", coef), Span::call_site());
                //                //quote! { #plus #coef_lit  #i_expr * #j_expr}
                //            } else if coef < 0 {
                //                panic!("Non-unity coefficient");
                //                //let coef_lit =
                //                //    LitInt::new(&format!("{}", -coef), Span::call_site());
                //                //quote! { #minus #coef_lit  self.#i_expr * r.#j_expr}
                //            } else {
                //                panic!("Bad coefficient (zero?)");
                //            }
                //        })
                //        .collect()
                //};

                let destructure = |oldvar, newvar, obj: &Object| {
                    if obj.is_scalar {
                        let newvar_ident = coefficient_symbol(newvar, &[]);
                        quote! { let #newvar_ident = #oldvar; }
                    } else {
                        let obj_type_name = obj.type_name_colons();
                        let obj_field_list: TokenStream = basis
                            .iter()
                            .zip(obj.select_components.iter())
                            .filter_map(|(b, is_selected)| {
                                is_selected.then(|| {
                                    let field = coefficient_symbol("a", b);
                                    let specific_newvar = coefficient_symbol(newvar, b);
                                    if field == specific_newvar {
                                        quote! { #field, }
                                    } else {
                                        quote! { #field: #specific_newvar, }
                                    }
                                })
                            })
                            .collect();

                        quote! { let #obj_type_name { #obj_field_list } = #oldvar; }
                    }
                };

                let destructure_self = destructure(quote! { self }, "a", lhs_obj);
                let destructure_rhs = destructure(quote! { r }, "b", rhs_obj);

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
                            let field = coefficient_symbol("a", b);
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
                    impl < T: ScalarRing > #op_trait < #rhs_type_name >  for #lhs_type_name {
                        #output_type

                        fn #op_fn (self, r: #rhs_type_name) -> #output_type_name {
                            #destructure_self
                            #destructure_rhs
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
                                let identifier = coefficient_symbol("a", b);
                                quote! {
                                    #identifier: T,
                                }
                            })
                        })
                        .collect();
                    let name = &obj.name;
                    quote! {
                        struct #name <T: ScalarRing> {
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
                false,
            );

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
                false,
            );

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
                false,
            );

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
                false,
            );

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
                false,
            );

            let op_trait = quote! { Transform };
            let op_fn = Ident::new("transform", Span::call_site());
            let sandwich_product = |i: usize, j: usize| {
                // Note: j is the bread of the sandwich
                let reverse_coef = |ix: usize| match (basis[ix].len() / 2) % 2 {
                    0 => 1,
                    1 => -1,
                    _ => panic!("Expected parity to be 0 or 1"),
                };

                // TODO this is wrong! it does not include cross terms
                let (coef_1, ix_int) = multiplication_table[j][i];
                let (coef_2, ix_result) = multiplication_table[ix_int][j];
                let coef = coef_1 * coef_2 * reverse_coef(j);
                assert!(ix_result == i);

                (coef, i)
            };
            /*let transform_code = gen_binary_operator(
                &basis,
                &objects,
                op_trait,
                op_fn,
                &obj,
                |a, b| generate_symbolic_product(a, b, sandwich_product),
                true,
            );*/
 // XXX

            quote! {
                // ===========================================================================
                // #name
                // ===========================================================================
                #wedge_product_code
                #vee_product_code
                #dot_product_code
                #geometric_product_code
                #commutator_product_code
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
