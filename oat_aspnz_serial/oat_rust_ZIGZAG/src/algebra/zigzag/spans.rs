use std::cell::Cell;
use std::collections::HashMap;
use std::hash::Hash;
use std::convert::TryInto;

use crate::algebra::chains::factored::{factor_boundary_matrix, FactoredBoundaryMatrix, JordanColumns};
use crate::algebra::chains::jordan;
use crate::algebra::matrices::operations::solve::triangle::{TriangularSolverMajorAscend, TriangularSolverMinorDescend};
use crate::algebra::rings::operator_structs::field_prime_order::BooleanFieldOperator;
use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;

    
use crate::algebra::rings::operator_traits::DivisionRing;
use crate::topology::simplicial::simplices::vector::subsimplices_dim_d_iter_descend;
use crate::topology::simplicial::{from::relation::BoundaryMatrixDowker, simplices::vector::{subsimplices_dim_0_thru_d_iter_ascend_dim_descend_lex, subsimplices_dim_d_iter_ascend}};
use crate::utilities::order::{is_sorted_strictly, OrderOperatorAuto, OrderOperatorByKey, OrderOperatorByKeyCustom};
use crate::utilities::sequences_and_ordinals::SortedVec;


use itertools::Itertools;
use serde::{Serialize, Deserialize};













/// Takes a collection of hyperedges and a ring operator as input; returns a factored boundary matrix as output.
pub fn factor_dowker_complex< RingOperator, Coefficient >( 
        dowker_simplices:           Vec<Vec<usize>>,
        ring_operator:              RingOperator,
        max_homology_dimension:     usize,    
    )  
    -> 
    FactoredBoundaryMatrix< 
        BoundaryMatrixDowker< usize, RingOperator, Coefficient >, 
        RingOperator, 
        OrderOperatorByKeyCustom< 
            Vec<usize>, 
            Coefficient, 
            (Vec<usize>,Coefficient), 
            OrderOperatorAuto,
        >, 
        Vec<usize>,
        (Vec<usize>, Coefficient),                    
        Vec<Vec<usize>>,
    > 

    where
        RingOperator:       Clone + DivisionRing< Coefficient >,
        Coefficient:        Clone + std::fmt::Debug + PartialEq,
{

    // Parameters
    // ----------
        
    // We will build a dowker complex.
    // A dowker complex is defined by a vertex set V and a family S
    // of subsets of V.  A subset of V forms a simplex iff it is 
    // a subset of some element of S.  We refer to the elements 
    // of S as "dowker simplices".
        
    // Each dowker simplex is represented by a SortedVec of vertices.
    // We store the list of all such simplices inside a larger vector.
    let dowker_simplices
        =   dowker_simplices
                .into_iter()
                .map( |x| SortedVec::new(x).unwrap() )  // we unwrap because `new` can return an error
                .collect_vec();
        
    //  Boundary matrix
    //  ---------------
        
    // This is a lazy object that generates rows/columns of the boundary matrix, on demand.
    let boundary_matrix = BoundaryMatrixDowker::new( dowker_simplices.clone(), ring_operator.clone() );
        
        
    // An iterator that runs over simplices of dimension 0 through max_homology_dimension,
    // ordered first by dimension (ascending) and second by lexicographic order (descending)
    let row_indices = boundary_matrix.row_indices_in_descending_order( max_homology_dimension as isize ).collect::<Vec<_>>();
        
    //  Homology computation (by matrix factorization)
    //  ----------------------------------------------
        
    // Return the factored boundary matrix
    factor_boundary_matrix( 
        boundary_matrix, 
        ring_operator, 
        OrderOperatorAuto, 
        row_indices,
    )
    
}








#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Span< Coefficient > {
    pub left_morphism:      VecOfVec< usize, Coefficient >,
    pub right_morphism:     VecOfVec< usize, Coefficient >,
    // pub left_grading:       Vec< usize >,
    // pub center_grading:     Vec< usize >,
    // pub right_grading:      Vec< usize >,
    // pub left_basis:         VecOfVec< Vec<usize>, Coefficient >,// Vec< ( Vec<usize>, Coefficient) >,
    pub left_basis_opt:     Option< VecOfVec< Vec<usize>, Coefficient > >,
    pub center_basis:       VecOfVec< Vec<usize>, Coefficient >,// Vec< ( Vec<usize>, Coefficient) >,
    pub right_basis:        VecOfVec< Vec<usize>, Coefficient >,// Vec< ( Vec<usize>, Coefficient) >,        
}








/// Returns matrices which represent linear operations on **ROWS**
pub fn induced_span< RingOperator, Coefficient >( 
        hypergraph_a:               Vec<Vec<usize>>, 
        hypergraph_b:               Vec<Vec<usize>>,
        ring_operator:              RingOperator,
        max_homology_dimension:     usize,
        return_left_basis:          bool,
    )
    ->
    Span< Coefficient >
    where
        RingOperator:       Clone + DivisionRing< Coefficient >,
        Coefficient:        Clone + std::fmt::Debug + PartialEq,


{

    let mut hypergraph_u                                =   hypergraph_a.clone();
    hypergraph_u.extend( hypergraph_b.iter().cloned() );


    // factor the boundary matrices
    let factored_boundary_matrix_a          =   factor_dowker_complex(hypergraph_a, ring_operator.clone(), max_homology_dimension );
    let factored_boundary_matrix_b          =   factor_dowker_complex(hypergraph_b, ring_operator.clone(), max_homology_dimension );
    let factored_boundary_matrix_u          =   factor_dowker_complex(hypergraph_u, ring_operator.clone(), max_homology_dimension );






    

    let homology_basis_in_u                             =   factored_boundary_matrix_u.basis_cohomology_vec_of_vec();



    let mut induced_maps        =   Vec::new();
    // let mut homology_cycle_bases                        =   Vec::new();

    // for each of the two bookend spaces (a and b), compute the induced map into u
    for factored_bookend_boundary_matrix in vec![ & factored_boundary_matrix_a, & factored_boundary_matrix_b ] {



        // place the homology basis vectors (or rather their indices) for space u into bijection with some integers
        //
        // NB: when OAT enumerates these indices it does an exhaustive search over the indices of the rows of the boundary matrix.
        //     this search visits simplices in the same order that we use when we factor the boundary matrix using the 
        //     cohomology algorithm. this means that we visit simplices in ascending order of dimension but *descending*
        //     lexicographic order within a dimension.
        //     therefore if simplices S and T have equal dimension and S ≤ T lexicographically, then  we will havee
        //     harmonic_basis_vector_index_in_bookend_to_ordinal( T ) ≤ harmonic_basis_vector_index_in_bookend_to_ordinal ( S )
        let mut harmonic_basis_vector_index_in_bookend_to_ordinal               =   HashMap::new();        
        
        for ( counter, index ) in factored_bookend_boundary_matrix.indices_harmonic().enumerate() {
            harmonic_basis_vector_index_in_bookend_to_ordinal.insert( index, counter  );
        }




        let cohomology_basis_in_u                       =     factored_boundary_matrix_u.basis_cohomology_vec_of_vec();


        // initialize the matrix that represents the induced map
        let mut induced_map                                 =   Vec::new();
    
        // iterate over basis vectors for the homology of the bookend space
        for cohomology_basis_vector_in_u        in    cohomology_basis_in_u.vec_of_vec() {

            let restriction_to_bookend                  =   cohomology_basis_vector_in_u
                                                                    .iter()
                                                                    .filter(
                                                                        |(simplex, coefficient)|
                                                                        {
                                                                            let simplex     =   SortedVec::new(simplex.clone() ).unwrap();
                                                                            let hypergraph  =   factored_bookend_boundary_matrix
                                                                                                    .umatch()
                                                                                                    .mapping_ref()
                                                                                                    .dowker_sets();
                                                                            for hyperedge in hypergraph {
                                                                                if hyperedge.contains_subset( &simplex ) { return true }
                                                                            }
                                                                            return false
                                                                        }
                                                                    );

            // find the family of coefficients a_i such that harmonic_basis_vector_in_bookend = sum_i ( a_i  * the_ith_jordan_basis_vector_in_u )
            let corresponding_linear_combination_in_bookend       =   TriangularSolverMajorAscend::solve(
                                                                    restriction_to_bookend.cloned(),
                                                                    factored_bookend_boundary_matrix.jordan_basis_matrix_dual(),
                                                                    ring_operator.clone(),
                                                                    OrderOperatorByKey::new(),
                                                                );

            // remove the coefficients that correspond to non-harmonic basis vectors (i.e. basis vectors in the coboundary)
            let mut projection                              =   corresponding_linear_combination_in_bookend
                                                                    .filter_map( 
                                                                        | (simplex,coefficient) |
                                                                        {   
                                                                            harmonic_basis_vector_index_in_bookend_to_ordinal
                                                                                .get( & simplex )
                                                                                .cloned()
                                                                                .map( | ordinal | (ordinal, coefficient) )
                                                                        }
    
                                                                    )
                                                                    .collect::<Vec<_>>(); 

            projection.reverse(); // reversing entries places them in sorted order                                                          
            

            // push the combination as a row to the induced map
            // 
            // NB: if we had not re-indexed the projection with integers then we sould have had to reverse the vector, because TriangularsolveMinorDescend
            //     returns entries in descending order of index. however, as we noted in the comments above, the map harmonic_basis_vector_index_in_bookend_to_ordinal
            //     reverses the order of indices within a fixed dimension
            induced_map.push( projection );

        }

        // convert the induced map matrix into a VecOfVec, and add it to our list
        let induced_map                                 =   VecOfVec::new(induced_map);
        induced_maps.push( induced_map );
    }


    let right_morphism                                  =   induced_maps.remove(1);
    let left_morphism                                   =   induced_maps.remove(0);   
    // let right_basis                                     =   homology_cycle_bases.remove(1);
    let left_basis_opt                                  =   if return_left_basis {  
                                                                                                    Some( factored_boundary_matrix_a.basis_cohomology_vec_of_vec() )
                                                                                            } else {
                                                                                                None
                                                                                            }; 


    Span{
        left_morphism,
        right_morphism,
        left_basis_opt,
        center_basis:       homology_basis_in_u,
        right_basis:        factored_boundary_matrix_b.basis_cohomology_vec_of_vec(),
    }

}






#[cfg(test)]
mod tests {    

    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

    #[test]
    fn doc_test_induced_cospan_two_cycles() {


        //
        // 3               5
        // 2       1       4
        //         0    
        //
        // the underlying simplicial complex we want to build looks like the result of
        // glueing two length-4 cycle graphs together by identifying a single edge
        // from each cycle

        // the "harmonic" basis of cycle representatives for the union (with vectors 
        // appearing in the same order in which they are returned by OAT iterators) is 
        // (omitting coefficients since we are working over the integers modulo 2):
        //
        // basis vector 0:  [0]
        // basis vector 1:  [0,1], [0,4], [1,5], [4,5]
        // basis vector 2:  [0,1], [0,2], [1,3], [2,3]        
        // 
        // this is because the iterator for harmonic indices follows the same order
        // as the iterator that runs over row indices; this iterator runs over indices
        // in ascending order of dimension (first) and descending lexicographic order
        // (second, for simplices of the same dimension)
        //
        // the "harmonic indicex" associated with each vector is the last index of
        // any nonzero entry in the vector. so they are:
        //
        // basis vector 0:  [0]
        // basis vector 1:  [4,5]
        // basis vector 2:  [2,3]
        //
        // as you can see, the indices in this sequence appear in sorted order
        // according to the rule outlined above
        //
        // the basis of cycle representatives for the other hypergraphs are
        //
        // hypergraph a:
        //   basis vector 0:    [0]
        //   basis vector 1:    [0,1], [0,2], [1,3], [2,3]
        // hypergraph b:
        //   basis vector 0:    [0]
        //   basis vector 1:    [0,1], [0,4], [1,5], [4,5]
        //
        // therefore the matrix representations for the maps induced by inclusion are
        //
        // hypergraph a:
        //   1 0 0
        //   0 0 1
        // hypergraph b:
        //   1 0 0
        //   0 1 0        

    


        // first half-filled cycle
        let hypergraph_a    =   vec![
                                    vec![ 0, 1 ],
                                    vec![ 0, 2 ],
                                    vec![ 1, 3 ],
                                    vec![ 2, 3 ],                                    
                                ];

        // second half-filled cycle                                
        let hypergraph_b    =   vec![
                                    vec![ 0, 1 ],
                                    vec![ 0, 4 ],
                                    vec![ 1, 5 ],
                                    vec![ 4, 5 ],
                                ];               

        // work over the two element field
        let ring_operator       =   BooleanFieldOperator::new();

        // set max homology dimension
        let max_homology_dimension                      =   2;

        // compute the span
        let return_left_basis               =   false;
        let span              =   induced_span( hypergraph_a, hypergraph_b, ring_operator, max_homology_dimension, return_left_basis );

        // check the matrices
        let matrix_a            =   & span.left_morphism;
        let matrix_b            =   & span.right_morphism;        
        let matrix_a_ground_truth: VecOfVec< usize, bool >   =   VecOfVec::new( 
                                            vec![
                                                vec![ (0,true)                        ], // this is for the 0-dimensional homology class
                                                vec![                       (2,true ) ], // this is for the "lefthand" cycle
                                            ]     
                                        );
        let matrix_b_ground_truth: VecOfVec< usize, bool >   =   VecOfVec::new( 
                                            vec![
                                                vec![ (0,true)                        ], // this is for the 0-dimensional homology class
                                                vec![           (1,true )             ], // this is for the "righthand" cycle
                                            ]     
                                        );                                        
        
        assert_eq!(  matrix_a, & matrix_a_ground_truth );
        assert_eq!(  matrix_b, & matrix_b_ground_truth );       
    }

}