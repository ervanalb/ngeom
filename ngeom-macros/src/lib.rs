extern crate proc_macro;
use core::cmp::Ordering;
use core::ops::{AddAssign, Mul};
use proc_macro2::Span;
use proc_macro2::TokenStream;
use quote::{quote, ToTokens, TokenStreamExt};
use syn::parse::{Error, Parse, ParseStream, Result};
use syn::punctuated::Punctuated;
use syn::{
    parse2, AttrStyle, Attribute, Fields, FieldsNamed, Ident, Item, ItemMacro, LitInt, Meta, Path,
    Token,
};

struct MultivectorStruct {
    ident: Ident,
    components: Vec<Ident>,
}

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

fn sign_from_parity(swaps: usize) -> isize {
    match swaps % 2 {
        0 => 1,
        1 => -1,
        _ => panic!("Expected parity to be 0 or 1"),
    }
}

#[derive(Default, Clone)]
struct SymbolicSumExpr(Vec<SymbolicProdExpr>);

#[derive(PartialEq, Eq, Clone)]
struct SymbolicProdExpr(isize, Vec<Symbol>);

#[derive(PartialOrd, Ord, PartialEq, Eq, Clone)]
enum Symbol {
    Scalar(Ident),
    StructField(Ident, Ident), // Two Idents: a var and a field (i.e. var.field)
}

impl ToTokens for Symbol {
    fn to_tokens(&self, tokens: &mut TokenStream) {
        match self {
            Symbol::StructField(var, field) => {
                var.to_tokens(tokens);
                <Token![.]>::default().to_tokens(tokens);
                field.to_tokens(tokens);
            }
            Symbol::Scalar(var) => {
                var.to_tokens(tokens);
            }
        }
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

impl Mul<isize> for SymbolicProdExpr {
    type Output = SymbolicProdExpr;
    fn mul(mut self, r: isize) -> SymbolicProdExpr {
        let SymbolicProdExpr(l_coef, _) = &mut self;
        *l_coef *= r;
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
            tokens.append_all(quote! { T::default() });
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
                        tokens.append_all(quote! { T::default() });
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
                        tokens.append_all(quote! { T::default() * });
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
        SymbolicSumExpr(l.into_iter().map(|lp| lp * r.clone()).collect())
    }
}

impl Mul<isize> for SymbolicSumExpr {
    type Output = SymbolicSumExpr;
    fn mul(self, r: isize) -> SymbolicSumExpr {
        let SymbolicSumExpr(l) = self;
        SymbolicSumExpr(l.into_iter().map(|lp| lp * r).collect())
    }
}

// Returns the right_complement of a basis element as a pair of
// (coef, complement_ix)
fn right_complement(right_complement_signs: &Vec<isize>, coef: isize, i: usize) -> (isize, usize) {
    let complement_ix = right_complement_signs.len() - i - 1;
    (coef * right_complement_signs[i], complement_ix)
}

// Returns the inverse of the right complement of a basis element
// (i.e. the left complement)
// as a pair of (coef, complement_ix)
fn left_complement(right_complement_signs: &Vec<isize>, coef: isize, i: usize) -> (isize, usize) {
    let complement_ix = right_complement_signs.len() - i - 1;
    (coef * right_complement_signs[complement_ix], complement_ix)
}

// We will represent a multivector as an array of coefficients on the basis elements.
// e.g. in 2D, there are 1 + 3 + 3 + 1 = 8 basis elements,
// and a full multivector uses all of them: [1, 1, 1, 1, 1, 1, 1, 1]
// An object such as a Bivector would only use a few of them: [0, 0, 0, 0, 1, 1, 1, 0]

#[derive(PartialEq)]
enum Object {
    Scalar,
    Struct(StructObject),
}

#[derive(PartialEq)]
struct StructObject {
    name: Ident,
    select_components: Vec<Option<(Ident, isize)>>,
    is_compound: bool,
}

impl Object {
    fn type_name(&self) -> TokenStream {
        match self {
            Object::Scalar => quote! { T },
            Object::Struct(StructObject { name, .. }) => {
                quote! { #name < T > }
            }
        }
    }
    fn type_name_colons(&self) -> TokenStream {
        match self {
            Object::Scalar => quote! { T },
            Object::Struct(StructObject { name, .. }) => {
                quote! { #name :: < T > }
            }
        }
    }
    fn has_component(&self, i: usize) -> bool {
        match self {
            Object::Scalar => i == 0,
            Object::Struct(StructObject {
                select_components, ..
            }) => select_components[i].is_some(),
        }
    }
    fn is_compound(&self) -> bool {
        match self {
            Object::Scalar => false,
            Object::Struct(StructObject { is_compound, .. }) => *is_compound,
        }
    }
    fn select_components(&self, var: Ident, len: usize) -> Vec<Option<(Symbol, isize)>> {
        match self {
            Object::Scalar => {
                let mut result = vec![None; len];
                result[0] = Some((Symbol::Scalar(var), 1));
                result
            }
            Object::Struct(StructObject {
                select_components, ..
            }) => select_components
                .iter()
                .map(|select_component| {
                    select_component.as_ref().map(|(field, coef)| {
                        (Symbol::StructField(var.clone(), field.clone()), *coef)
                    })
                })
                .collect(),
        }
    }
}

// Function to generate a simple unary operator on an object,
// where entries can be rearranged and coefficients can be modified
// (e.g. for implementing neg, right_complement, reverse)
// The result will be a list of symbolic expressions,
// one for each coefficient value in the resultant multivector.
fn generate_symbolic_rearrangement<F: Fn(isize, usize) -> (isize, usize)>(
    select_components: &[Option<(Symbol, isize)>],
    op: F,
) -> Vec<SymbolicSumExpr> {
    // Generate the sum
    let mut result: Vec<SymbolicSumExpr> = vec![Default::default(); select_components.len()];

    for (i, is_selected) in select_components.iter().enumerate() {
        if let Some((symbol, coef)) = is_selected {
            let (coef, result_basis_ix) = op(*coef, i);
            result[result_basis_ix] += SymbolicProdExpr(coef, vec![symbol.clone()]);
        }
    }

    result.into_iter().map(|expr| expr.simplify()).collect()
}

fn generate_symbolic_norm<F: Fn(isize, usize, isize, usize) -> (isize, usize)>(
    select_components: &[Option<(Symbol, isize)>],
    product: F,
    sqrt: bool,
) -> Vec<SymbolicSumExpr> {
    // Generate the product
    let mut expressions: Vec<SymbolicSumExpr> = vec![Default::default(); select_components.len()];
    for (i, (i_symbol, i_coef)) in select_components
        .iter()
        .enumerate()
        .filter_map(|(i, selected)| selected.as_ref().map(|selected| (i, selected)))
    {
        for (j, (j_symbol, j_coef)) in select_components
            .iter()
            .enumerate()
            .filter_map(|(j, selected)| selected.as_ref().map(|selected| (j, selected)))
        {
            let (coef, ix) = product(*i_coef, i, *j_coef, j);

            expressions[ix] += SymbolicProdExpr(coef, vec![i_symbol.clone(), j_symbol.clone()]);
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

// Function to generate a sum of two objects, e.g. for overloading + or -
// The result will be a list of symbolic expressions,
// one for each coefficient value in the resultant multivector.
fn generate_symbolic_sum(
    select_components_a: &[Option<(Symbol, isize)>],
    select_components_b: &[Option<(Symbol, isize)>],
    coef_a: isize,
    coef_b: isize,
) -> Vec<SymbolicSumExpr> {
    // Generate the sum
    let mut result: Vec<SymbolicSumExpr> = vec![Default::default(); select_components_a.len()];

    for (i, (a_selected, b_selected)) in select_components_a
        .iter()
        .zip(select_components_b.iter())
        .enumerate()
    {
        if let Some((symbol_a, coef_symbol_a)) = a_selected {
            result[i] += SymbolicProdExpr(coef_symbol_a * coef_a, vec![symbol_a.clone()]);
        }
        if let Some((symbol_b, coef_symbol_b)) = b_selected {
            result[i] += SymbolicProdExpr(coef_symbol_b * coef_b, vec![symbol_b.clone()]);
        }
    }

    result.into_iter().map(|expr| expr.simplify()).collect()
}

// Function to generate a product of two objects, e.g. geometric product, wedge product, etc.
// The result will be a list of symbolic expressions,
// one for each coefficient value in the resultant multivector.
fn generate_symbolic_product<F: Fn(isize, usize, isize, usize) -> (isize, usize)>(
    select_components_a: &[Option<(Symbol, isize)>],
    select_components_b: &[Option<(Symbol, isize)>],
    product: F,
) -> Vec<SymbolicSumExpr> {
    // Generate the product
    let mut result: Vec<SymbolicSumExpr> = vec![Default::default(); select_components_a.len()];
    for (i, (i_symbol, i_coef)) in select_components_a
        .iter()
        .enumerate()
        .filter_map(|(i, selected)| selected.as_ref().map(|selected| (i, selected)))
    {
        for (j, (j_symbol, j_coef)) in select_components_b
            .iter()
            .enumerate()
            .filter_map(|(j, selected)| selected.as_ref().map(|selected| (j, selected)))
        {
            let (coef, ix) = product(*i_coef, i, *j_coef, j);

            result[ix] += SymbolicProdExpr(coef, vec![i_symbol.clone(), j_symbol.clone()]);
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
    F1: Fn(isize, usize, isize, usize) -> (isize, usize),
    F2: Fn(isize, usize, isize, usize) -> (isize, usize),
>(
    select_components_a: &[Option<(Symbol, isize)>],
    select_components_b: &[Option<(Symbol, isize)>],
    product_1: F1,
    product_2: F2,
) -> Vec<SymbolicSumExpr> {
    // Generate the first intermediate product B PRODUCT1 A
    // where i maps to components of B, and j maps to components of A
    let mut intermediate_result: Vec<SymbolicSumExpr> =
        vec![Default::default(); select_components_a.len()];
    for (i, (i_symbol, i_coef)) in select_components_b
        .iter()
        .enumerate()
        .filter_map(|(i, selected)| selected.as_ref().map(|selected| (i, selected)))
    {
        for (j, (j_symbol, j_coef)) in select_components_a
            .iter()
            .enumerate()
            .filter_map(|(j, selected)| selected.as_ref().map(|selected| (j, selected)))
        {
            let (coef, ix) = product_1(*i_coef, i, *j_coef, j);
            intermediate_result[ix] +=
                SymbolicProdExpr(coef, vec![i_symbol.clone(), j_symbol.clone()]);
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
        for (j, (j_symbol, j_coef)) in select_components_b
            .iter()
            .enumerate()
            .filter_map(|(j, selected)| selected.as_ref().map(|selected| (j, selected)))
        {
            let (coef, ix) = product_2(1, i, *j_coef, j);
            let new_term = SymbolicProdExpr(coef, vec![j_symbol.clone()]);
            let result_term = intermediate_term.clone() * new_term;
            result[ix] += result_term;
        }
    }

    result.into_iter().map(|expr| expr.simplify()).collect()
}

fn find_output_object<'a>(
    objects: &'a [Object],
    output_expressions: &[SymbolicSumExpr],
) -> Option<&'a Object> {
    let select_output_components: Vec<_> = output_expressions
        .iter()
        .map(|SymbolicSumExpr(e)| e.len() != 0)
        .collect();
    objects.iter().find(|o| {
        select_output_components
            .iter()
            .enumerate()
            .find(|&(i, &out_c)| out_c && !o.has_component(i))
            .is_none()
    })
}

fn gen_unary_operator(
    objects: &[Object],
    op_trait: TokenStream,
    op_fn: Ident,
    obj: &Object,
    expressions: &[SymbolicSumExpr],
    alias: Option<(TokenStream, Ident)>,
) -> TokenStream {
    if matches!(obj, Object::Scalar) {
        // Do not generate operations with the scalar being the LHS--
        // typically because these would violate rust's orphan rule
        // or result in conflicting trait implementations
        return quote! {};
    };

    // Figure out what the type of the output is
    let output_object = find_output_object(&objects, &expressions);

    let Some(output_object) = output_object else {
        // No output object matches the result we got,
        // so don't generate any code
        return quote! {};
    };

    if matches!(output_object, Object::Scalar) && expressions[0].0.len() == 0 {
        // This operation unconditionally returns 0,
        // so invoking it is probably a type error--
        // do not generate code for it
        return quote! {};
    }

    let output_type_name = &output_object.type_name_colons();
    let type_name = &obj.type_name();

    let return_expr = match output_object {
        Object::Scalar => {
            let expr = &expressions[0];
            quote! { #expr }
        }
        Object::Struct(output_struct_object) => {
            let output_fields: TokenStream = output_struct_object
                .select_components
                .iter()
                .zip(expressions.iter())
                .map(|(select_component, expr)| {
                    if let Some((field, coef)) = select_component {
                        let expr = expr.clone() * *coef;
                        quote! { #field: #expr, }
                    } else {
                        return quote! {};
                    }
                })
                .collect();
            quote! {
                #output_type_name {
                    #output_fields
                }
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

fn gen_binary_operator<
    F: Fn(&[Option<(Symbol, isize)>], &[Option<(Symbol, isize)>]) -> Vec<SymbolicSumExpr>,
>(
    basis_element_count: usize,
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
            if matches!(lhs_obj, Object::Scalar) {
                // Do not generate operations with the scalar being the LHS--
                // these violate rust's orphan rule
                // or result in conflicting trait implementations.
                // Technically we could allow this for custom ops such as T::cross(r: Vector<T>)
                // but we won't for the sake of consistency.
                // Use Vector<T>::cross(r: T) instead.
                return quote! {};
            };

            let expressions = op(
                &lhs_obj
                    .select_components(Ident::new("self", Span::call_site()), basis_element_count),
                &rhs_obj.select_components(Ident::new("r", Span::call_site()), basis_element_count),
            );

            // Figure out what the type of the output is
            let output_object = find_output_object(&objects, &expressions);

            let Some(output_object) = output_object else {
                // No output object matches the result we got,
                // so don't generate any code
                return quote! {};
            };

            let rhs_type_name = &rhs_obj.type_name();
            let lhs_type_name = &lhs_obj.type_name();
            let output_type_name = &output_object.type_name_colons();

            if !implicit_promotion_to_compound
                && output_object.is_compound()
                && !(lhs_obj.is_compound() || rhs_obj.is_compound())
            {
                // Do not create compound objects unintentionally.
                // Only allow returning a compound object
                // when taking products of compound objects and other objects
                return quote! {};
            }

            if matches!(output_object, Object::Scalar) && expressions[0].0.len() == 0 {
                // This operation unconditionally returns 0,
                // so invoking it is probably a type error--
                // do not generate code for it.

                return quote! {};
            }

            let return_expr = match output_object {
                Object::Scalar => {
                    let expr = &expressions[0];
                    quote! { #expr }
                }
                Object::Struct(output_struct_object) => {
                    let output_fields: TokenStream = output_struct_object
                        .select_components
                        .iter()
                        .zip(expressions.iter())
                        .map(|(select_component, expr)| {
                            if let Some((field, coef)) = select_component {
                                let expr = expr.clone() * *coef;
                                quote! { #field: #expr, }
                            } else {
                                return quote! {};
                            }
                        })
                        .collect();
                    quote! {
                        #output_type_name {
                            #output_fields
                        }
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

fn implement_geometric_algebra(
    basis_vector_idents: Vec<Ident>,
    metric: Vec<isize>,
    multivector_structs: Vec<MultivectorStruct>,
) -> Result<TokenStream> {
    // Sanity checks
    if basis_vector_idents.len() == 0 {
        return Err(Error::new(Span::call_site(), "Basis vector set is empty"));
    }
    if basis_vector_idents.len() != metric.len() {
        return Err(Error::new(
            Span::call_site(),
            "Metric and basis are different sizes",
        ));
    }
    if multivector_structs.len() == 0 {
        return Err(Error::new(
            Span::call_site(),
            "No multivector structs defined",
        ));
    }

    // The number of dimensions in the algebra
    // (e.g. use a 4D algebra to represent 3D geometry)
    let dimension = metric.len();

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
        basis
    };

    let basis_element_count = basis.len();

    // Analyze the basis vector identifiers to find the common prefix,
    // as well as the part that varies
    let (ident_prefix, ident_variants) = {
        let first = basis_vector_idents
            .first()
            .map(|x| x.to_string())
            .unwrap_or_default();

        let len = first.chars().count();
        assert!(len >= 1, "Identifier should be non-zero length");
        let prefix = first.chars().take(len - 1).collect::<String>();

        let variants = basis_vector_idents
            .iter()
            .map(|basis_ident| {
                let basis_ident_str = basis_ident.to_string();
                if basis_ident_str.chars().count() != len || !basis_ident_str.starts_with(&prefix) {
                    return Err(Error::new(
                        basis_ident.span(),
                        "Bad identifier name: must be common prefix + 1 char",
                    ));
                }

                Ok(basis_ident_str.chars().nth(len - 1).unwrap())
            })
            .collect::<Result<Vec<_>>>()?;
        (prefix, variants)
    };
    let ident_prefix_len = ident_prefix.chars().count();

    // Start with the scalar type
    let mut objects = vec![Object::Scalar];
    let struct_objects = multivector_structs
        .iter()
        .map(|multivector_struct| {
            let name = multivector_struct.ident.clone();
            let mut select_components = vec![None; basis_element_count];

            // Put together the select_components array
            for component_ident in multivector_struct.components.iter() {
                let component_ident_str = component_ident.to_string();

                if component_ident_str.chars().count() <= ident_prefix_len
                    || !component_ident_str.starts_with(&ident_prefix)
                {
                    return Err(Error::new(
                        component_ident.span(),
                        "Bad identifier: must be common prefix + 1 or more chars",
                    ));
                }
                let variant_product = component_ident_str
                    .chars()
                    .skip(ident_prefix_len)
                    .collect::<Vec<_>>();
                let mut b = variant_product
                    .iter()
                    .map(|c| {
                        ident_variants
                            .iter()
                            .position(|c2| c2 == c)
                            .ok_or(Error::new(
                                component_ident.span(),
                                "Bad identifier: must be composed of basis vectors",
                            ))
                    })
                    .collect::<Result<Vec<_>>>();

                // Don't immediately throw the error if we get something that is not composed of basis vectors--we will special-case the scalar value
                // which can be any single unused letter if there is no common prefix.
                if b.is_err()
                    && select_components[0].is_none()
                    && ident_prefix_len == 0
                    && component_ident_str.chars().count() == 1
                {
                    b = Ok(vec![]);
                }

                let mut b = b?;

                let sign = sign_from_parity(bubble_sort_count_swaps(&mut b));
                let basis_index = basis.iter().position(|b2| b2 == &b).ok_or(Error::new(
                    component_ident.span(),
                    "Bad identifier: cannot repeat basis vectors",
                ))?;

                if select_components[basis_index].is_some() {
                    return Err(Error::new(
                        component_ident.span(),
                        "Bad identifier: duplicate",
                    ));
                }
                select_components[basis_index] = Some((component_ident.clone(), sign));
            }

            // Check if object is mixed-grade
            let grades = select_components
                .iter()
                .enumerate()
                .filter_map(|(i, component)| component.as_ref().map(|_| basis[i].len()))
                .collect::<Vec<_>>();
            let is_compound = grades.windows(2).any(|g| g[0] != g[1]);

            Ok(Object::Struct(StructObject {
                name,
                select_components,
                is_compound,
            }))
        })
        .collect::<Result<Vec<_>>>()?;

    objects.extend(struct_objects);

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
            sign_from_parity(bubble_sort_count_swaps(product.as_mut()))
        })
        .collect();

    // Generate dot product multiplication table for the basis elements.
    // The dot product always produces a scalar.
    let dot_product_multiplication_table: Vec<Vec<isize>> = {
        let multiply_basis_vectors = |ei: usize, ej: usize| {
            // This would need to get more complicated for CGA
            // where the metric g is not a diagonal matrix
            if ei == ej {
                metric[ei]
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

    let dot_product_f = |coef_i: isize, i: usize, coef_j: isize, j: usize| {
        let coef_mul = dot_product_multiplication_table[i][j];
        (coef_i * coef_j * coef_mul, 0)
    };

    let anti_dot_product_f = |coef_i: isize, i: usize, coef_j: isize, j: usize| {
        let (i_coef, i) = right_complement(&right_complement_signs, coef_i, i);
        let (j_coef, j) = right_complement(&right_complement_signs, coef_j, j);
        let (coef, ix) = dot_product_f(i_coef, i, j_coef, j);
        let (coef, ix) = left_complement(&right_complement_signs, coef, ix);
        (coef, ix)
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
                    coef *= metric[e];
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

    let geometric_product_f = |coef_i: isize, i: usize, coef_j: isize, j: usize| {
        let (coef_mul, ix) = geometric_product_multiplication_table[i][j];
        (coef_i * coef_j * coef_mul, ix)
    };

    let geometric_antiproduct_f = |coef_i: isize, i: usize, coef_j: isize, j: usize| {
        let (coef_i, i) = right_complement(&right_complement_signs, coef_i, i);
        let (coef_j, j) = right_complement(&right_complement_signs, coef_j, j);
        let (coef, ix) = geometric_product_f(coef_i, i, coef_j, j);
        let (coef, ix) = left_complement(&right_complement_signs, coef, ix);
        (coef, ix)
    };

    // TODO Fix handling of antiscalar!
    let anti_scalar_ops: TokenStream = {
        let field = Ident::new(
            &(ident_prefix + &ident_variants.into_iter().collect::<String>()),
            Span::call_site(),
        );

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
        }
    };

    let impl_code: TokenStream = objects
        .iter()
        .map(|obj| {
            // Derive From / Into
            let from_code: TokenStream = {
                if matches!(obj, Object::Scalar) {
                    // Do not implement From on the scalar--
                    // typically because these would violate rust's orphan rule
                    // and are also not useful
                    return quote! {};
                }

                objects.iter().map(|other_obj| {
                    // See if the object we are converting from
                    // contains a non-empty strict subset of the components in our object

                    let is_subset = (0..basis_element_count).all(|i| obj.has_component(i) || !other_obj.has_component(i));
                    if !is_subset { return quote! {}; }

                    let is_same_object = (0..basis_element_count).all(|i| obj.has_component(i) == other_obj.has_component(i));
                    if is_same_object { return quote! {}; }

                    let is_not_empty = (0..basis_element_count).any(|i| obj.has_component(i) && other_obj.has_component(i));
                    if !is_not_empty { return quote! {}; }

                    let my_type_name = obj.type_name();
                    let my_type_name_colons = obj.type_name_colons();
                    let other_type_name = other_obj.type_name();

                    let expressions: Vec<_> = other_obj.select_components(Ident::new("value", Span::call_site()), basis_element_count).iter().map(|select_component| {
                        match select_component {
                            Some((symbol, coef)) => SymbolicSumExpr(vec![SymbolicProdExpr(*coef, vec![symbol.clone()])]),
                            None => Default::default(),
                        }
                    }).collect();

                    let return_expr = match obj {
                        Object::Scalar => {
                            let expr = &expressions[0];
                            quote! { #expr }
                        }
                        Object::Struct(output_struct_object) => {
                            let output_fields: TokenStream = output_struct_object
                                .select_components
                                .iter()
                                .zip(expressions.iter())
                                .map(|(select_component, expr)| {
                                    if let Some((field, coef)) = select_component {
                                        let expr = expr.clone() * *coef;
                                        quote! { #field: #expr, }
                                    } else {
                                        return quote! {};
                                    }
                                })
                                .collect();
                            quote! {
                                #my_type_name_colons {
                                    #output_fields
                                }
                            }
                        }
                    };

                    quote! {
                        impl<T: core::default::Default> From<#other_type_name> for #my_type_name {
                            fn from(value: #other_type_name) -> #my_type_name {
                                #return_expr
                            }
                        }
                    }
                }).collect()
            };

            let obj_self_components = &obj.select_components(Ident::new("self", Span::call_site()), basis_element_count);

            // Overload unary -
            let op_trait = quote! { core::ops::Neg };
            let op_fn = Ident::new("neg", Span::call_site());
            let neg_code = gen_unary_operator(
                &objects,
                op_trait,
                op_fn,
                &obj,
                &generate_symbolic_rearrangement(&obj_self_components, |coef: isize, i: usize| (-coef, i)),
                None,
            );

            // Add a method A.reverse()
            let op_trait = quote! { Reverse };
            let op_fn = Ident::new("reverse", Span::call_site());
            let reverse_f = |coef: isize, i: usize| {
                    let coef_rev = sign_from_parity((basis[i].len() / 2) % 2);
                    (coef * coef_rev, i)
                };
            let reverse_expressions = generate_symbolic_rearrangement(&obj_self_components, reverse_f);
            let reverse_code = gen_unary_operator(&objects, op_trait, op_fn, &obj, &reverse_expressions, None);

            // Add a method A.anti_reverse()
            let op_trait = quote! { AntiReverse };
            let op_fn = Ident::new("anti_reverse", Span::call_site());
            let anti_reverse_f = |coef: isize, i: usize| {
                    let (coef, i) = right_complement(&right_complement_signs, coef, i);
                    let (coef, i) = reverse_f(coef, i);
                    let (coef, i) = left_complement(&right_complement_signs, coef, i);
                    (coef, i)
                };
            let anti_reverse_expressions = generate_symbolic_rearrangement(&obj_self_components, anti_reverse_f);
            let anti_reverse_code = gen_unary_operator(&objects, op_trait, op_fn, &obj, &anti_reverse_expressions, None);

            // Add a method A.right_complement()
            let op_trait = quote! { RightComplement };
            let op_fn = Ident::new("right_complement", Span::call_site());
            let right_complement_expressions = generate_symbolic_rearrangement(&obj_self_components, |coef: isize, i: usize| right_complement(&right_complement_signs, coef, i));
            let right_complement_code = gen_unary_operator(&objects, op_trait, op_fn, &obj, &right_complement_expressions, None);

            // Add a method A.left_complement()
            let op_trait = quote! { LeftComplement };
            let op_fn = Ident::new("left_complement", Span::call_site());
            let left_complement_expressions = generate_symbolic_rearrangement(&obj_self_components, |coef: isize, i: usize| left_complement(&right_complement_signs, coef, i));
            let left_complement_code = gen_unary_operator(&objects, op_trait, op_fn, &obj, &left_complement_expressions, None);

            // Add a method A.bulk()
            let op_trait = quote! { Bulk };
            let op_fn = Ident::new("bulk", Span::call_site());
            let bulk_f = |coef: isize, i: usize| {
                // This would need to be changed for CGA where the dot product multiplication table
                // is not diagonal
                let (zero_or_one, _) = dot_product_f(1, i, 1, i);
                (coef * zero_or_one, i)
            };
            let bulk_expressions = generate_symbolic_rearrangement(&obj_self_components, bulk_f);
            let bulk_code = gen_unary_operator(&objects, op_trait, op_fn, &obj, &bulk_expressions, None);

            // Add a method A.weight()
            let op_trait = quote! { Weight };
            let op_fn = Ident::new("weight", Span::call_site());
            let weight_f = |coef: isize, i: usize| {
                let (coef, i) = right_complement(&right_complement_signs, coef, i);
                let (coef, i) = bulk_f(coef, i);
                let (coef, i) = left_complement(&right_complement_signs, coef, i);
                (coef, i)
            };
            let weight_expressions = generate_symbolic_rearrangement(&obj_self_components, weight_f);
            let weight_code = gen_unary_operator(&objects, op_trait, op_fn, &obj, &weight_expressions, None);

            // Add a method A.bulk_dual() which computes A★
            let op_trait = quote! { BulkDual };
            let op_fn = Ident::new("bulk_dual", Span::call_site());
            let bulk_dual_f = |coef: isize, i: usize| {
                let (coef, i) = bulk_f(coef, i);
                let (coef, i) = right_complement(&right_complement_signs, coef, i);
                (coef, i)
            };
            let bulk_dual_expressions = generate_symbolic_rearrangement(&obj_self_components, bulk_dual_f);
            let bulk_dual_code =
                gen_unary_operator(&objects, op_trait, op_fn, &obj, &bulk_dual_expressions, None);

            // Add a method A.weight_dual() which computes A☆
            let op_trait = quote! { WeightDual };
            let op_fn = Ident::new("weight_dual", Span::call_site());
            let alias_trait = quote! { Normal };
            let alias_fn = Ident::new("normal", Span::call_site());
            let weight_dual_f = |coef: isize, i: usize| {
                let (coef, i) = weight_f(coef, i);
                let (coef, i) = right_complement(&right_complement_signs, coef, i);
                (coef, i)
            };
            let weight_dual_expressions = generate_symbolic_rearrangement(&obj_self_components, weight_dual_f);
            let weight_dual_code =
                gen_unary_operator(&objects, op_trait, op_fn, &obj, &weight_dual_expressions, Some((alias_trait, alias_fn)));


            // Add a method A.bulk_norm_squared()
            let op_trait = quote! { BulkNormSquared };
            let op_fn = Ident::new("bulk_norm_squared", Span::call_site());
            // Squared norm uses the product A • A

            let bulk_norm_squared_expressions = generate_symbolic_norm(&obj_self_components, dot_product_f, false);
            let bulk_norm_squared_code =
                gen_unary_operator(&objects, op_trait, op_fn, &obj, &bulk_norm_squared_expressions, None);

            let bulk_norm_code = if !bulk_norm_squared_code.is_empty() {
                // Add a method A.bulk_norm() if it is possible to symbolically take the sqrt
                let op_trait = quote! { BulkNorm };
                let op_fn = Ident::new("bulk_norm", Span::call_site());
                let bulk_norm_expressions = generate_symbolic_norm(&obj_self_components, dot_product_f, true);
                let bulk_norm_code = gen_unary_operator(&objects, op_trait, op_fn, &obj, &bulk_norm_expressions, None);
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
            let weight_norm_squared_expressions = generate_symbolic_norm(&obj_self_components, anti_dot_product_f, false);

            let weight_norm_squared_code =
                gen_unary_operator(&objects, op_trait, op_fn, &obj, &weight_norm_squared_expressions, None);

            let weight_norm_code = if !weight_norm_squared_code.is_empty() {
                // Add a method A.weight_norm() if it is possible to symbolically take the sqrt
                let op_trait = quote! { WeightNorm };
                let op_fn = Ident::new("weight_norm", Span::call_site());
                let weight_norm_expressions = generate_symbolic_norm(&obj_self_components, anti_dot_product_f, true);
                let weight_norm_code =
                    gen_unary_operator(&objects, op_trait, op_fn, &obj, &weight_norm_expressions, None);
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
            let hat_code = if !matches!(obj, Object::Scalar) {
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
                        Self: WeightNorm<Output=AntiScalar<T>> // TODO: Figure out anti-scalar
                                                               // handling
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
                basis_element_count,
                &objects,
                op_trait,
                op_fn,
                &obj,
                |a, b| generate_symbolic_sum(a, b, 1, 1),
                true,  // implicit_promotion_to_compound
                None, // alias
            );

            // Overload -
            let op_trait = quote! { core::ops::Sub };
            let op_fn = Ident::new("sub", Span::call_site());
            let sub_code = gen_binary_operator(
                basis_element_count,
                &objects,
                op_trait,
                op_fn,
                &obj,
                |a, b| generate_symbolic_sum(a, b, 1, -1),
                true,  // implicit_promotion_to_compound
                None, // alias
            );

            // Add a method A.wedge(B) which computes A ∧ B
            let op_trait = quote! { Wedge };
            let op_fn = Ident::new("wedge", Span::call_site());
            let alias_trait = quote! { Join };
            let alias_fn = Ident::new("join", Span::call_site());
            let wedge_product_f = |coef_i: isize, i: usize, coef_j: isize, j: usize| {
                let (coef, ix) = geometric_product_f(coef_i, i, coef_j, j);
                // Select grade s + t
                let s = basis[i].len();
                let t = basis[j].len();
                let u = basis[ix].len();
                let coef = if s + t == u { coef } else { 0 };
                (coef, ix)
            };
            let wedge_product_code = gen_binary_operator(
                basis_element_count,
                &objects,
                op_trait,
                op_fn,
                &obj,
                |a, b| generate_symbolic_product(a, b, wedge_product_f),
                false, // implicit_promotion_to_compound
                Some((alias_trait, alias_fn)), // alias
            );

            // Add a method A.anti_wedge(B) which computes A ∨ B
            let op_trait = quote! { AntiWedge };
            let op_fn = Ident::new("anti_wedge", Span::call_site());
            let alias_trait = quote! { Meet };
            let alias_fn = Ident::new("meet", Span::call_site());
            let anti_wedge_product_f = |coef_i: isize, i: usize, coef_j: isize, j: usize| {
                let (coef_i, i) = right_complement(&right_complement_signs, coef_i, i);
                let (coef_j, j) = right_complement(&right_complement_signs, coef_j, j);

                let (coef, ix) = geometric_product_f(coef_i, i, coef_j, j);
                // Select grade s + t
                let s = basis[i].len();
                let t = basis[j].len();
                let u = basis[ix].len();
                let coef = if s + t == u { coef } else { 0 };

                let (coef, ix) = left_complement(&right_complement_signs, coef, ix);
                (coef, ix)
            };
            let anti_wedge_product_code = gen_binary_operator(
                basis_element_count,
                &objects,
                op_trait,
                op_fn,
                &obj,
                |a, b| generate_symbolic_product(a, b, anti_wedge_product_f),
                false, // implicit_promotion_to_compound
                Some((alias_trait, alias_fn)), // alias
            );

            // Add a method A.dot(B) which computes A • B
            let op_trait = quote! { Dot };
            let op_fn = Ident::new("dot", Span::call_site());

            let dot_product_code = gen_binary_operator(
                basis_element_count,
                &objects,
                op_trait,
                op_fn,
                &obj,
                |a, b| generate_symbolic_product(a, b, dot_product_f),
                false, // implicit_promotion_to_compound
                None, // alias
            );

            // Add a method A.anti_dot(B) which computes A ∘ B
            let op_trait = quote! { AntiDot };
            let op_fn = Ident::new("anti_dot", Span::call_site());

            let anti_dot_product_code = gen_binary_operator(
                basis_element_count,
                &objects,
                op_trait,
                op_fn,
                &obj,
                |a, b| generate_symbolic_product(a, b, anti_dot_product_f),
                false, // implicit_promotion_to_compound
                None, // alias
            );

            // Implement the geometric product ⟑
            let op_trait = quote! { WedgeDot };
            let op_fn = Ident::new("wedge_dot", Span::call_site());

            let wedge_dot_product_code = gen_binary_operator(
                basis_element_count,
                &objects,
                op_trait,
                op_fn,
                &obj,
                |a, b| generate_symbolic_product(a, b, geometric_product_f),
                false, // implicit_promotion_to_compound
                None, // alias
            );

            // Implement the geometric antiproduct ⟇
            let op_trait = quote! { AntiWedgeDot };
            let op_fn = Ident::new("anti_wedge_dot", Span::call_site());
            let alias_trait = quote! { Compose };
            let alias_fn = Ident::new("compose", Span::call_site());

            let anti_wedge_dot_product_code = gen_binary_operator(
                basis_element_count,
                &objects,
                op_trait,
                op_fn,
                &obj,
                |a, b| generate_symbolic_product(a, b, geometric_antiproduct_f),
                false, // implicit_promotion_to_compound
                Some((alias_trait, alias_fn)), // alias
            );

            // Overload * for scalar multiplication
            let op_trait = quote! { core::ops::Mul };
            let op_fn = Ident::new("mul", Span::call_site());
            let scalar_product = |coef_i: isize, i: usize, coef_j: isize, j: usize| {
                let (coef, ix) = geometric_product_f(coef_i, i, coef_j, j);
                // Require s or t be a scalar
                let s = basis[i].len();
                let t = basis[j].len();
                let coef = if s == 0 || t == 0 { coef } else { 0 };
                (coef, ix)
            };
            let scalar_product_code = gen_binary_operator(
                basis_element_count,
                &objects,
                op_trait,
                op_fn,
                &obj,
                |a, b| generate_symbolic_product(a, b, scalar_product),
                true, // implicit_promotion_to_compound
                None, // alias
            );

            // Implement A.anti_mul(B) for anti-scalar multiplication
            let op_trait = quote! { AntiMul };
            let op_fn = Ident::new("anti_mul", Span::call_site());
            let anti_scalar_f = |coef_i: isize, i: usize, coef_j: isize, j: usize|  {
                let (coef_i, i) = right_complement(&right_complement_signs, coef_i, i);
                let (coef_j, j) = right_complement(&right_complement_signs, coef_j, j);

                let (coef, ix) = scalar_product(coef_i, i, coef_j, j);

                let (coef, ix) = left_complement(&right_complement_signs, coef, ix);
                (coef, ix)
            };
            let anti_scalar_product_code = gen_binary_operator(
                basis_element_count,
                &objects,
                op_trait,
                op_fn,
                &obj,
                |a, b| generate_symbolic_product(a, b, anti_scalar_f),
                true, // implicit_promotion_to_compound
                None, // alias
            );

            // Add a method A.bulk_expansion(B) which computes A ∧ B★
            let op_trait = quote! { BulkExpansion };
            let op_fn = Ident::new("bulk_expansion", Span::call_site());
            let bulk_expansion_f = |coef_i: isize, i: usize, coef_j: isize, j: usize| {
                let (coef_j, j) = bulk_dual_f(coef_j, j);
                let (coef, ix) = wedge_product_f(coef_i, i, coef_j, j);
                (coef, ix)
            };
            let bulk_expansion_code = gen_binary_operator(
                basis_element_count,
                &objects,
                op_trait,
                op_fn,
                &obj,
                |a, b| generate_symbolic_product(a, b, bulk_expansion_f),
                false, // implicit_promotion_to_compound
                None, // alias
            );

            // Add a method A.weight_expansion(B) which computes A ∧ B☆
            let op_trait = quote! { WeightExpansion };
            let op_fn = Ident::new("weight_expansion", Span::call_site());
            let alias_trait = quote! { SupersetOrthogonalTo };
            let alias_fn = Ident::new("superset_orthogonal_to", Span::call_site());
            let weight_expansion_f = |coef_i: isize, i: usize, coef_j: isize, j: usize| {
                let (coef_j, j) = weight_dual_f(coef_j, j);
                let (coef, ix) = wedge_product_f(coef_i, i, coef_j, j);
                (coef, ix)
            };
            let weight_expansion_code = gen_binary_operator(
                basis_element_count,
                &objects,
                op_trait,
                op_fn,
                &obj,
                |a, b| generate_symbolic_product(a, b, weight_expansion_f),
                false, // implicit_promotion_to_compound
                Some((alias_trait, alias_fn)), // alias
            );

            // Add a method A.bulk_contraction(B) which computes A ∨ B★
            let op_trait = quote! { BulkContraction };
            let op_fn = Ident::new("bulk_contraction", Span::call_site());
            let bulk_contraction_f = |coef_i: isize, i: usize, coef_j: isize, j: usize| {
                let (coef_j, j) = bulk_dual_f(coef_j, j);
                let (coef, ix) = anti_wedge_product_f(coef_i, i, coef_j, j);
                (coef, ix)
            };
            let bulk_contraction_code = gen_binary_operator(
                basis_element_count,
                &objects,
                op_trait,
                op_fn,
                &obj,
                |a, b| generate_symbolic_product(a, b, bulk_contraction_f),
                false, // implicit_promotion_to_compound
                None, // alias
            );

            // Add a method A.weight_contraction(B) which computes A ∨ B☆
            let op_trait = quote! { WeightContraction };
            let op_fn = Ident::new("weight_contraction", Span::call_site());
            let alias_trait = quote! { SubsetOrthogonalTo };
            let alias_fn = Ident::new("subset_orthogonal_to", Span::call_site());
            let weight_contraction_f = |coef_i: isize, i: usize, coef_j: isize, j: usize| {
                let (coef_j, j) = weight_dual_f(coef_j, j);
                let (coef, ix) = anti_wedge_product_f(coef_i, i, coef_j, j);
                (coef, ix)
            };
            let weight_contraction_code = gen_binary_operator(
                basis_element_count,
                &objects,
                op_trait,
                op_fn,
                &obj,
                |a, b| generate_symbolic_product(a, b, weight_contraction_f),
                false, // implicit_promotion_to_compound
                Some((alias_trait, alias_fn)), // alias
            );

            // Add a method A.anti_commutator(B) which computes (A ⟇ B - B ⟇ A) / 2
            let op_trait = quote! { AntiCommutator };
            let op_fn = Ident::new("anti_commutator", Span::call_site());
            let commutator_product = |coef_i: isize, i: usize, coef_j: isize, j: usize| {
                let (coef_1, ix) = geometric_antiproduct_f(coef_i, i, coef_j, j);
                let (coef_2, ix_2) = geometric_antiproduct_f(coef_j, j, coef_i, i);

                let coef = coef_1 - coef_2;
                assert!(ix == ix_2);
                assert!(coef % 2 == 0);
                let coef = coef / 2;
                (coef, ix)
            };
            let commutator_product_code = gen_binary_operator(
                basis_element_count,
                &objects,
                op_trait,
                op_fn,
                &obj,
                |a, b| generate_symbolic_product(a, b, commutator_product),
                false, // implicit_promotion_to_compound
                None, // alias
            );

            // Add a method A.transform(B) which computes B̰ ⟇ A ⟇ B
            let op_trait = quote! { Transform };
            let op_fn = Ident::new("transform", Span::call_site());
            let transform_1 = |coef_i: isize, i: usize, coef_j: isize, j: usize| {
                // Compute first half of B̰ ⟇ A ⟇ B
                // Where i maps to B, and j maps to A.
                // In part 2, we will compute the geometric antiproduct of this intermediate result
                // with B

                let (coef_i, i) = anti_reverse_f(coef_i, i);
                geometric_antiproduct_f(coef_i, i, coef_j, j)
            };
            // Compute second half of B̰ ⟇ A ⟇ B
            // In part 1, we computed the intermediate result B̰ ⟇ A which maps to i here.
            // j maps to B.
            let transform_2 = geometric_antiproduct_f;

            let transform_code = gen_binary_operator(
                basis_element_count,
                &objects,
                op_trait,
                op_fn,
                &obj,
                |a, b| {
                    generate_symbolic_double_product(
                        a,
                        b,
                        transform_1,
                        transform_2,
                    )
                },
                false, // implicit_promotion_to_compound
                None, // alias
            );

            // Add a method A.reverse_transform(B) which computes B ⟇ A ⟇ B̰
            let op_trait = quote! { ReverseTransform };
            let op_fn = Ident::new("reverse_transform", Span::call_site());
            // Compute first half of B ⟇ A ⟇ B̰
            // Where i maps to B, and j maps to A.
            // In part 2, we will compute the geometric antiproduct of this intermediate result
            // with B̰
            let reverse_transform_1 = geometric_antiproduct_f;
            let reverse_transform_2 = |coef_i: isize, i: usize, coef_j: isize, j: usize| {
                // Compute second half of B ⟇ A ⟇ B̰
                // In part 1, we computed the intermediate result B ⟇ A which maps to i here.
                // j maps to B.

                let (coef_j, j) = anti_reverse_f(coef_j, j);
                geometric_antiproduct_f(coef_i, i, coef_j, j)
            };
            let reverse_transform_code = gen_binary_operator(
                basis_element_count,
                &objects,
                op_trait,
                op_fn,
                &obj,
                |a, b| {
                    generate_symbolic_double_product(
                        a,
                        b,
                        reverse_transform_1,
                        reverse_transform_2,
                    )
                },
                false, // implicit_promotion_to_compound
                None, // alias
            );

            // Implement A.projection(B) which computes B ∨ (A ∧  B☆)
            let op_trait = quote! { Projection };
            let op_fn = Ident::new("projection", Span::call_site());
            let projection_product_1 = |coef_i: isize, i: usize, coef_j: isize, j: usize| {
                // Compute second half of B ∨ (A ∧ B☆)
                // Where i maps to B, and j maps to A.
                // In part 2, we will compute the geometric product of B with this intermediate result
                weight_expansion_f(coef_j, j, coef_i, i)
            };
            let projection_product_2 = |coef_i: isize, i: usize, coef_j: isize, j: usize| {
                // Compute second half of B ∨ (A ∧ B☆)
                // In part 1, we computed the intermediate result A ∧ B☆ which maps to i here.
                // j maps to B.
                anti_wedge_product_f(coef_j, j, coef_i, i)
            };
            let projection_code = gen_binary_operator(
                basis_element_count,
                &objects,
                op_trait,
                op_fn,
                &obj,
                |a, b| {
                    generate_symbolic_double_product(
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
            let anti_projection_product_1 = |coef_i: isize, i: usize, coef_j: isize, j: usize| {
                // Compute second half of B ∧ (A ∨ B☆)
                // Where i maps to B, and j maps to A.
                // In part 2, we will compute the geometric product of B with this intermediate result
                weight_contraction_f(coef_j, j, coef_i, i)
            };
            let anti_projection_product_2 = |coef_i: isize, i: usize, coef_j: isize, j: usize| {
                // Compute second half of B ∧ (A ∨ B☆)
                // In part 1, we computed the intermediate result A ∨ B☆ which maps to i here.
                // j maps to B.
                wedge_product_f(coef_j, j, coef_i, i)
            };
            let anti_projection_code = gen_binary_operator(
                basis_element_count,
                &objects,
                op_trait,
                op_fn,
                &obj,
                |a, b| {
                    generate_symbolic_double_product(
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
            let central_projection_product_1 = |coef_i: isize, i: usize, coef_j: isize, j: usize| {
                // Compute second half of B ∨ (A ∧ B★)
                // Where i maps to B, and j maps to A.
                // In part 2, we will compute the geometric product of B with this intermediate result
                bulk_expansion_f(coef_j, j, coef_i, i)
            };
            let central_projection_product_2 = |coef_i: isize, i: usize, coef_j: isize, j: usize| {
                // Compute second half of B ∨ (A ∧ B★)
                // In part 1, we computed the intermediate result A ∧ B★ which maps to i here.
                // j maps to B.
                anti_wedge_product_f(coef_j, j, coef_i, i)
            };
            let central_projection_code = gen_binary_operator(
                basis_element_count,
                &objects,
                op_trait,
                op_fn,
                &obj,
                |a, b| {
                    generate_symbolic_double_product(
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
            let central_anti_projection_product_1 = |coef_i: isize, i: usize, coef_j: isize, j: usize| {
                // Compute second half of B ∧ (A ∨ B★)
                // Where i maps to B, and j maps to A.
                // In part 2, we will compute the geometric product of B with this intermediate result
                bulk_contraction_f(coef_j, j, coef_i, i)
            };
            let central_anti_projection_product_2 = |coef_i: isize, i: usize, coef_j: isize, j: usize| {
                // Compute second half of B ∧ (A ∨ B★)
                // In part 1, we computed the intermediate result A ∨ B★ which maps to i here.
                // j maps to B.
                wedge_product_f(coef_j, j, coef_i, i)
            };
            let central_anti_projection_code = gen_binary_operator(
                basis_element_count,
                &objects,
                op_trait,
                op_fn,
                &obj,
                |a, b| {
                    generate_symbolic_double_product(
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
            let motor_to_f = |coef_i: isize, i: usize, coef_j: isize, j: usize| {
                let (coef_i, i) = anti_reverse_f(coef_i, i);
                geometric_antiproduct_f(coef_i, i, coef_j, j)
            };

            let motor_to_code = gen_binary_operator(
                basis_element_count,
                &objects,
                op_trait,
                op_fn,
                &obj,
                |a, b| generate_symbolic_product(a, b, motor_to_f),
                true, // implicit_promotion_to_compound
                None, // alias
            );

            quote! {
                // ===========================================================================
                // #name
                // ===========================================================================

                #from_code
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
                #reverse_transform_code
                #motor_to_code
            }
        })
        .collect();

    Ok(quote! {
        #anti_scalar_ops
        #impl_code
    })
}

/// A wrapper around Vec<syn::Item> so we can implement Parse on it
struct VecItem(Vec<Item>);

impl Parse for VecItem {
    fn parse(input: ParseStream) -> Result<Self> {
        let mut items = Vec::<Item>::new();
        while !input.is_empty() {
            items.push(input.parse()?);
        }
        Ok(VecItem(items))
    }
}

struct BasisVectorIdents(Vec<Ident>);

impl Parse for BasisVectorIdents {
    fn parse(input: ParseStream) -> Result<Self> {
        let ident_list = Punctuated::<Ident, Token![,]>::parse_terminated(input)?;
        Ok(BasisVectorIdents(ident_list.into_iter().collect()))
    }
}

struct Metric(Vec<isize>);

impl Parse for Metric {
    fn parse(input: ParseStream) -> Result<Self> {
        let lit_list = Punctuated::<LitInt, Token![,]>::parse_terminated(input)?;
        Ok(Metric(
            lit_list
                .into_iter()
                .map(|lit| lit.base10_parse::<isize>())
                .collect::<Result<Vec<_>>>()?,
        ))
    }
}

fn geometric_algebra2(code: TokenStream) -> Result<TokenStream> {
    let VecItem(mut items) = parse2(code)?;

    let mut metric = None;
    let mut basis = None;
    let mut multivector_structs = Vec::<MultivectorStruct>::new();

    // Idents to match when parsing
    let multivector_ident = Ident::new("multivector", Span::call_site());
    let metric_ident = Ident::new("metric", Span::call_site());
    let basis_ident = Ident::new("basis", Span::call_site());

    let mut err: Result<()> = Ok(());

    let mut append_err = |new_e: Error| {
        err = match &mut err {
            Ok(()) => Err(new_e),
            Err(old_e) => {
                old_e.combine(new_e);
                Err(old_e.clone())
            }
        };
    };

    items.retain_mut(|item| {
        match item {
            Item::Macro(ItemMacro {
                mac: item_macro, ..
            }) => {
                let macro_ident = item_macro.path.get_ident();
                if macro_ident == Some(&basis_ident) {
                    if basis.is_some() {
                        append_err(Error::new(
                            macro_ident.unwrap().span(),
                            "Duplicate basis definition",
                        ));
                    } else {
                        // Parse the basis as a bracket-enclosed comma-separated list of identifiers
                        let parsed_basis: Result<BasisVectorIdents> = item_macro.parse_body();
                        match parsed_basis {
                            Ok(BasisVectorIdents(parsed_basis)) => {
                                basis = Some(parsed_basis);
                            }
                            Err(e) => {
                                append_err(e);
                            }
                        }
                    }
                    false // Do not retain the basis! macro
                } else if macro_ident == Some(&metric_ident) {
                    if metric.is_some() {
                        append_err(Error::new(
                            macro_ident.unwrap().span(),
                            "Duplicate metric definition",
                        ));
                    } else {
                        // Parse the metric as a bracket-enclosed comma-separated list of identifiers
                        let parsed_metric: Result<Metric> = item_macro.parse_body();
                        match parsed_metric {
                            Ok(Metric(parsed_metric)) => {
                                metric = Some(parsed_metric);
                            }
                            Err(e) => {
                                append_err(e);
                            }
                        }
                    }
                    false // Do not retain the metric! macro
                } else {
                    true // Retain unrecognized macros
                }
            }
            Item::Struct(item_struct) => {
                let mut has_multivector_attribute = false;

                item_struct.attrs.retain(|attr| {
                    if let Attribute {
                        style: AttrStyle::Outer,
                        meta: Meta::Path(Path { segments, .. }),
                        ..
                    } = attr
                    {
                        // Do not retain the #[multivector] attribute
                        // but make note that we found it
                        if segments.len() == 1
                            && segments.first().unwrap().ident == multivector_ident
                        {
                            has_multivector_attribute = true;
                            false
                        } else {
                            true
                        }
                    } else {
                        true
                    }
                });

                if has_multivector_attribute {
                    // This struct is a multivector!
                    // Parse its fields and add it to the struct_info list
                    if let Fields::Named(FieldsNamed { named: fields, .. }) = &item_struct.fields {
                        // TODO check types
                        match fields
                            .iter()
                            .map(|field| Ok(field.ident.as_ref().unwrap().clone()))
                            .collect::<Result<Vec<_>>>()
                        {
                            Ok(components) => {
                                multivector_structs.push(MultivectorStruct {
                                    ident: item_struct.ident.clone(),
                                    components,
                                });
                            }
                            Err(e) => {
                                append_err(e);
                            }
                        }
                    } else {
                        append_err(Error::new(
                            item_struct.ident.span(),
                            "Multivector must have named fields",
                        ));
                    }
                }

                true // Retain structs
            }
            _ => {
                true // Retain unrecognized items
            }
        }
    });
    err?;

    let Some(basis) = basis else {
        return Err(Error::new(
            Span::call_site(),
            "Missing basis![..] definition",
        ));
    };

    let Some(metric) = metric else {
        return Err(Error::new(
            Span::call_site(),
            "Missing metric![..] definition",
        ));
    };

    // Generate code
    let generated_code = implement_geometric_algebra(basis, metric, multivector_structs)?;

    // Add generated code to original code
    let mut code = TokenStream::new();
    code.append_all(items);
    generated_code.to_tokens(&mut code);
    Ok(code)
}

#[proc_macro]
pub fn geometric_algebra(code: proc_macro::TokenStream) -> proc_macro::TokenStream {
    // Parse input
    geometric_algebra2(code.into())
        .unwrap_or_else(Error::into_compile_error)
        .into()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn derive_geometric_algebra() {
        let _result = geometric_algebra2(quote! {
            basis![w, x, y];
            metric![0, 1, 1];
            #[multivector]
            struct Vector<T> {
                x: T,
                y: T,
                w: T,
            }
            #[multivector]
            struct Bivector<T> {
                wx: T,
                wy: T,
                xy: T,
            }
            #[multivector]
            struct AntiScalar<T> {
                wxy: T,
            }
            #[multivector]
            struct AntiEven<T> {
                a: T,
                wx: T,
                wy: T,
                xy: T,
            }
            #[multivector]
            struct AntiOdd<T> {
                x: T,
                y: T,
                w: T,
                wxy: T,
            }
            #[multivector]
            struct Multivector<T> {
                a: T,
                x: T,
                y: T,
                w: T,
                wx: T,
                wy: T,
                xy: T,
                wxy: T,
            }
        })
        .unwrap();
        //println!("{}", result);
        //panic!();
    }
}
