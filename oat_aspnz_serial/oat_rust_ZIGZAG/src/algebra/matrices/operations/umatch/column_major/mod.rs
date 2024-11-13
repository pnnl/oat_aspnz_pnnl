



use itertools::Itertools;
use crate::algebra::matrices::operations::Umatch;


use crate::algebra::matrices::operations::solve::triangle::{TriangularSolverMinorDescend, TriangularSolverMajorAscend};
use crate::algebra::matrices::types::bimajor::{MatrixBimajor, MatrixBimajorData};
use crate::algebra::matrices::types::packet::MatrixAlgebraPacket;
use crate::algebra::matrices::types::transpose::AntiTranspose;
use crate::algebra::vectors::entries::{KeyValSet, KeyValNew};
use crate::algebra::matrices::operations::transform_entry_wise::ReindexSquareMatrix;
use crate::algebra::matrices::types::matching::{GeneralizedMatchingArrayWithMajorOrdinals};

use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
use crate::algebra::matrices::query::{ViewRowAscend, ViewColDescend, IndicesAndCoefficients, ViewsMinorDescend};
use crate::utilities::functions::evaluate::{ EvaluateFunctionFnMutWrapper };


use crate::algebra::matrices::operations::umatch::row_major::comb::{CombDomain, CombCodomain, CombDomainInv, CombCodomainInv};

use crate::utilities::iterators::general::{FilterOutMembers};
use crate::utilities::iterators::merge::hit::{HitMerge};
use crate::algebra::vectors::entries::{KeyValGet};
use crate::algebra::rings::operator_traits::{Semiring, Ring, DivisionRing};

use crate::utilities::order::{JudgePartialOrder, OrderOperatorByKeyCustom, ReverseOrder};


use crate::algebra::vectors::operations::{Scale, VectorOperations, Simplify};

use std::collections::HashMap;
use std::hash::Hash;
use std::fmt::Debug;
use std::iter::{Cloned, Peekable};



use crate::algebra::matrices::operations::solve::echelon::{EchelonSolverMajorAscendWithMajorKeys};

use crate::algebra::matrices::operations::transform_vector_wise::{OnlyKeyMinInsideCollection, OnlyKeyMinOutsideCollection, OnlyKeyMajInsideCollection, OnlyKeyMajOutsideCollection};




//  ===========================================================================






pub struct UmatchColumnMajor< Mapping, RingOperator, ReverseOrderOperatorRowEntries, ReverseOrderOperatorColEntries > 
    where   
        Mapping:                    ViewColDescend + IndicesAndCoefficients,
        Mapping::ColIndex:          Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::RowIndex:          Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::ViewMinorDescend:   IntoIterator,
        Mapping::EntryMinor:        KeyValGet< Mapping::RowIndex, Mapping::Coefficient >,
{
    row_major_umatch_of_antitranspose:                  Umatch<
                                                            AntiTranspose< Mapping >, 
                                                            RingOperator, 
                                                            ReverseOrderOperatorColEntries, 
                                                            ReverseOrderOperatorRowEntries, 
                                                        >
}








//  ---------------------------------------------------------------------------------------------------------
//  U-MATCH -- CONSTRUCTORS
//  ---------------------------------------------------------------------------------------------------------



impl < Mapping, RingOperator, OrderOperatorRowIndex, OrderOperatorColIndex, >  

    UmatchColumnMajor
        < 
            Mapping, 
            RingOperator, 
            OrderOperatorByKeyCustom< Mapping::ColIndex, Mapping::Coefficient, Mapping::EntryMajor, ReverseOrder< OrderOperatorColIndex >, >,  // recall that entries in the major views have minor keys
            OrderOperatorByKeyCustom< Mapping::RowIndex, Mapping::Coefficient, Mapping::EntryMinor, ReverseOrder< OrderOperatorRowIndex >, >, // recall that entries in the minor views have major keys
        >  

    where   
        Mapping::ColIndex:                  Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::RowIndex:                  Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct 
        Mapping:                            ViewColDescend + IndicesAndCoefficients,
        Mapping::ViewMinorDescend:          Clone + IntoIterator,
        Mapping::ViewMinorDescendIntoIter:  Clone,        
        Mapping::EntryMinor:                KeyValSet< Mapping::RowIndex, Mapping::Coefficient > + Debug + Clone,    // !!!!!!!!!!!!!! REMOVE THE DEBUG REQUIREMENT AFTER DEBUGGING!!!!!            
        Mapping::ColIndex:                 Clone + Hash + std::cmp::Eq + Debug, 
        Mapping::RowIndex:                 Clone + Hash + std::cmp::Eq + Debug,     // !!! remove Debug eventually       
        Mapping::Coefficient:                 Clone + Debug,            // !!! remove Debug eventually       
        RingOperator:           Clone + Semiring< Mapping::Coefficient > + Ring< Mapping::Coefficient > + DivisionRing< Mapping::Coefficient >,
        OrderOperatorColIndex:   Clone + JudgePartialOrder <  Mapping::ColIndex >, // !!! remove clone
        OrderOperatorRowIndex:   Clone + JudgePartialOrder <  Mapping::RowIndex >,
        HitMerge<Peekable<Scale<<Mapping::ViewMinorDescend as IntoIterator>::IntoIter, Mapping::RowIndex, RingOperator, Mapping::Coefficient>>, OrderOperatorByKeyCustom< Mapping::RowIndex, Mapping::Coefficient, Mapping::EntryMinor, OrderOperatorColIndex>>: Clone, // !!!! remove this                        

                
{

    /// Generate a new U-match factorization
    /// 
    /// # Arguments
    /// 
    /// - `mapping`: matrix you want to factor (or "mapping matrix")
    /// - `iter_keymin`: an iterator that runs over the column indices (= minor indices) of the matrix, in *strictly ascending order*
    /// - `ring_operator`: operator for the coefficient ring
    /// - `order_operator_keymaj`: order operator for the row indices (= major keys)
    /// - `order_operator_keymin`: order operator for the column indices (= minor keys)
    /// 
    /// # Special features of this function
    /// 
    /// Although the Umatch struct can take arbitrary orders the entries in its major and minor views, this
    /// constructor requires the user to provide an order on minor and major keys.  This is a safety measure, to prevent
    /// common forms of errors.  However, if you really want a custom order on entries, there are other ways to construct
    /// a U-match factorization, besides this constructor.
    pub fn factor
            < IterRowIndex > ( 
                mapping:                    Mapping, 
                iter_keymin:                IterRowIndex,
                ring_operator:              RingOperator,
                order_operator_keymaj:      OrderOperatorRowIndex,                   
                order_operator_keymin:      OrderOperatorColIndex,             
            ) 
        -> 
        Self

        where
            IterRowIndex:               Iterator < Item = Mapping::ColIndex >,        
    {

        let antitransposed_mapping                      =   AntiTranspose::new( mapping );

        let row_major_umatch_of_antitranspose           =   Umatch::factor(
                                                                antitransposed_mapping,
                                                                iter_keymin,  // argument name: iter_keymaj
                                                                ring_operator,
                                                                ReverseOrder::new( order_operator_keymin ), // argument name: order_operator_keymaj
                                                                ReverseOrder::new( order_operator_keymaj ), // argument name: order_operator_keymin
                                                            );

        UmatchColumnMajor{ row_major_umatch_of_antitranspose }
    }

}


impl < Mapping, RingOperator, ReverseOrderOperatorRowEntries, ReverseOrderOperatorColEntries, >  

    UmatchColumnMajor
    < Mapping, RingOperator, ReverseOrderOperatorRowEntries, ReverseOrderOperatorColEntries, >  
    
    where   
        Mapping::ColIndex:                  Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::RowIndex:                  Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct 
        Mapping:                            ViewColDescend + IndicesAndCoefficients,
        Mapping::ViewMinorDescend:           IntoIterator,
        Mapping::EntryMinor:                KeyValSet< Mapping::RowIndex, Mapping::Coefficient > + Debug + Clone,    // !!!!!!!!!!!!!! REMOVE THE DEBUG REQUIREMENT AFTER DEBUGGING!!!!!            
        Mapping::ColIndex:                 Clone + Hash + std::cmp::Eq + Debug, 
        Mapping::RowIndex:                 Clone + Hash + std::cmp::Eq + Debug,     // !!! remove Debug eventually       
        Mapping::Coefficient:                 Clone + Debug,            // !!! remove Debug eventually       
        RingOperator:           Clone + Semiring< Mapping::Coefficient > + Ring< Mapping::Coefficient > + DivisionRing< Mapping::Coefficient >,
{


    /// Returns a references to the row-major Umatch factorization of the anti-transpose of the boundary matrix
    pub fn row_major_umatch_of_antitranspose_ref( &self ) ->
        & Umatch<
            AntiTranspose< Mapping >, 
            RingOperator, 
            ReverseOrderOperatorColEntries, 
            ReverseOrderOperatorRowEntries, 
        >
    {
        & self.row_major_umatch_of_antitranspose
    }


    /// Returns a copy of the ring operator
    pub fn ring_operator( &self ) -> RingOperator 
        where
            RingOperator:   Clone,
    { self.row_major_umatch_of_antitranspose.ring_operator() }

    /// Returns a copy of the order comparator for index-value pairs whose index is a `Mapping::ColIndex`.
    pub fn order_operator_major( &self ) -> ReverseOrder< ReverseOrderOperatorRowEntries >
        where
            ReverseOrderOperatorRowEntries:   Clone,
    {  ReverseOrder::new( self.row_major_umatch_of_antitranspose.order_operator_minor() ) } 
    
    /// Returns a copy of the **inverted** order comparator for index-value pairs whose index is a `Mapping::ColIndex`.
    pub fn order_operator_major_reverse( &self ) -> ReverseOrderOperatorRowEntries
        where
            ReverseOrderOperatorRowEntries:   Clone,
    { self.row_major_umatch_of_antitranspose.order_operator_minor() }        

    /// Returns a copy of the order comparator for index-value pairs whose index is a `Mapping::RowIndex`.
    pub fn order_operator_minor( &self ) -> ReverseOrder< ReverseOrderOperatorColEntries >
        where
            ReverseOrderOperatorColEntries:   Clone,
    { ReverseOrder::new( self.row_major_umatch_of_antitranspose.order_operator_major() ) } 
    
    /// Returns a copy of the **inverted** order comparator for index-value pairs whose index is a `Mapping::RowIndex`.
    pub fn order_operator_minor_reverse( &self ) -> ReverseOrderOperatorColEntries
        where
            ReverseOrderOperatorColEntries:   Clone,
    { self.row_major_umatch_of_antitranspose.order_operator_major() }  


   

    /// Returns the (row-major) codomain COMB
    /// 
    /// # Design notes
    /// 
    /// This matrix cannot be indexed by integers without extra work, since we only assign ordinals to *matched* rows of the mapping array.
    pub fn comb_codomain( &self ) -> 
        AntiTranspose<
            CombDomainInv< '_, AntiTranspose<Mapping>, RingOperator, ReverseOrderOperatorColEntries, ReverseOrderOperatorRowEntries >
        >  
    {
        AntiTranspose::new(
            self.row_major_umatch_of_antitranspose.comb_domain_inv() 
        )
    }

    /// Returns the (row-major) inverse of the codomain COMB
    /// 
    /// # Design notes
    /// 
    /// This matrix cannot be indexed by integers without extra work, since we only assign ordinals to *matched* rows of the mapping array.
    pub fn comb_codomain_inv( &self ) -> 
        AntiTranspose<
            CombDomain< '_, AntiTranspose< Mapping >, RingOperator, ReverseOrderOperatorColEntries, ReverseOrderOperatorRowEntries >  
        >
    {
        AntiTranspose::new(        
            self.row_major_umatch_of_antitranspose.comb_domain()
        )
    }  
    
    /// Returns the (row-major) domain COMB, indexed by `Mapping::ColIndex`
    /// 
    /// # Design notes
    /// 
    /// This matrix cannot be indexed by integers without extra work, since we only assign ordinals to *matched* columns of the mapping array.
    pub fn comb_domain( &self ) -> 
        AntiTranspose<
            CombCodomainInv< '_, AntiTranspose< Mapping >, RingOperator, ReverseOrderOperatorColEntries, ReverseOrderOperatorRowEntries >
        >
    {
        AntiTranspose::new(
            self.row_major_umatch_of_antitranspose.comb_codomain_inv()
        )
    }

    /// Returns the (row-major) inverse of the domain COMB
    /// 
    /// # Design notes
    /// 
    /// This matrix cannot be indexed by integers without extra work, since we only assign ordinals to *matched* columns of the mapping array.
    pub fn comb_domain_inv( &self ) -> 
        AntiTranspose<
            CombCodomain< '_, AntiTranspose< Mapping >, RingOperator, ReverseOrderOperatorColEntries, ReverseOrderOperatorRowEntries >  
        >
    {
        AntiTranspose::new(
            self.row_major_umatch_of_antitranspose.comb_codomain()
        )
    }      

}    

