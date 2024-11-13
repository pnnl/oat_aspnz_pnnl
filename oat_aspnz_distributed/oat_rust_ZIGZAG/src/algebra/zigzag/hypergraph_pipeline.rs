use std::collections::HashSet;
use std::time::Instant;
use std::sync::Arc;

use lamellar::active_messaging::prelude::*;
use lamellar::darc::prelude::*;
use lamellar::{Serialize,Deserialize,de::DeserializeOwned,bincode};

use once_cell::sync::Lazy;

use itertools::Itertools;
use num::Integer;

use crate::algebra::chains::factored;
use crate::algebra::matrices::operations::multiply::vector_matrix_multiply_major_ascend_simplified;
use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
use crate::algebra::rings::operator_structs::field_prime_order::BooleanFieldOperator;
use crate::algebra::zigzag::cospan_pipeline;
use crate::algebra::zigzag::cospans::async_induced_cospan;
use crate::algebra::{
    chains::factored::factor_boundary_matrix, rings::operator_traits::DivisionRing,
};
use crate::algebra::zigzag::cospans::Cospan;
use crate::utilities::order::OrderOperatorByKey;

use super::decompose::Diagonalization;
use super::{
    cospans::{factor_dowker_complex, induced_cospan},
    decompose::{QuiverReprsentation, SingleBarBasisVectorIndexLedger},
};

use derive_getters::{Dissolve, Getters};
use derive_new::new;

// THIS JUST MAKES THE CODE EASIER TO READ
type Simplex = Vec<usize>;
type Chain<RingElement> = Vec<(Simplex, RingElement)>;

/// The dimension of a linear combination of simplices
///
/// Assumes that the linear combination is homogeous with respect to simplex dimension
pub fn chain_dimension<RingElement>(chain: &Vec<(Simplex, RingElement)>) -> usize {
    chain[0].0.len() - 1
}

#[derive(new, Getters, Clone, Debug, Dissolve, Eq, PartialEq, Ord, PartialOrd)]
pub struct IntervalSubmoduleBasis<RingElement> {
    pub leftmost_vertex: usize,
    pub basis_vectors: Vec<Chain<RingElement>>,
}

impl<RingElement> IntervalSubmoduleBasis<RingElement> {
    pub fn interval_endpoints(&self) -> (usize, usize) {
        let left = self.leftmost_vertex.clone();
        let right = left + self.basis_vectors.len();
        return (left, right);
    }
}

static WORLD: Lazy<lamellar::LamellarWorld> =
    Lazy::new(|| lamellar::LamellarWorldBuilder::new().build());

/// Returns an interval decomposition of the zigzag module of a sequence of hypergraphs connected by their unions.
///
/// Concretely, we take a sequence of hypergraphs `A, B, C, ..` and decompose the zigzag module for
/// `A --> A u B <-- B --> B u C <-- ..`.
///
/// The coefficient field is determined by the choice of `ring_operator`.
///
/// # Examples
///
/// ```
/// use oat_rust::algebra::rings::operator_structs::field_prime_order::BooleanFieldOperator;
/// use oat_rust::topology::simplicial::from::relation::sideways_ladder_edges;
/// use oat_rust::algebra::zigzag::hypergraph_pipeline::interval_decomposition_for_zigzag_of_hypgeraph_unions;
///
///
/// // Construct length-2 sequence of hypergraphs of the following form. (in this
/// // case there are actually just simple undirected graphs)
/// //
/// // hypergraph 0:
/// //
/// // 1 ---- 3 ---- 5
/// // |      |      |
/// // |      |      |
/// // 0 ---- 2 ---- 4
/// //
/// // hypergraph 1:
/// //
/// //        3 ---- 5 ---- 7
/// //        |      |      |
/// //        |      |      |
/// //        2 ---- 4 ---- 6
/// //
/// // --------------------------------------------------------------------
///
/// let number_of_ladders                           =   2;
/// let holes_per_ladder                            =   2;
///
/// let mut hypergraphs                             =   Vec::new();
/// for offset_from_left in 0 .. number_of_ladders {
///     hypergraphs.push(
///         sideways_ladder_edges(offset_from_left, holes_per_ladder)  
///     );   
/// }
///
/// // Compute the barcode for the corresponding zigzag:
/// //
/// // hypergraph_0 ----> ( hypergraph_0  U  hypergraph_1 ) <---- hypergraph_1
/// //
/// // --------------------------------------------------------------------
///
/// let ring_operator                               =   BooleanFieldOperator::new();
/// let max_homology_dimension                      =   2;
/// let interval_decomposition                      =   interval_decomposition_for_zigzag_of_hypgeraph_unions(hypergraphs, ring_operator, max_homology_dimension);
/// let barcode                                     =   interval_decomposition.iter()
///                                                         .map(| (dim, submodule)| ( dim.clone(), submodule.interval_endpoints() ) )
///                                                         .collect::<Vec<_>>();
///
/// // Verify that the barcode is correct
/// // Note that each pair (a,b) represents a half-open interval of form [a,b)
/// // The actual barcode is
/// //
/// //             
/// // dimension 0:
/// //       
/// //      0 ---- 1 ---- 2             
/// //
/// // dimension 1  
/// //      
/// //      0 ---- 1
/// //      0 ---- 1 ---- 2
/// //             1 ---- 2
/// //
/// // --------------------------------------------------------------------
///         
/// let barcode_ground_truth: Vec< (usize,(usize,usize)) >                        =   vec![
///                                                         ( 0, (0,3) ),
///                                                         ( 1, (0,2) ),
///                                                         ( 1, (0,3) ),
///                                                         ( 1, (1,3) ),                                                                                                                                                                                               
///                                                     ];
///
/// assert_eq!( barcode, barcode_ground_truth );
/// ```

#[lamellar::AmLocalData]
struct InducedCospanAm< RingOperator,RingElement> {
    hypergraphs: Arc<Vec<Vec<Vec<usize>>>>,
    a: usize,
    b: usize,
    ring_operator: RingOperator,
    max_homology_dimension: usize,
    _phantom: std::marker::PhantomData<RingElement>,
}

#[lamellar::local_am]
impl<RingOperator, RingElement> LamellarAM for InducedCospanAm< RingOperator, RingElement>
where
    RingOperator: Clone + Send + Sync  + DivisionRing<RingElement> + std::fmt::Debug + 'static,
    RingElement: Clone + Send + Sync + Serialize + DeserializeOwned + std::fmt::Debug + Eq + Ord + 'static, 
{
    async fn exec(self) -> Vec<u8>{ //Cospan<RingElement> {
        let cospan = induced_cospan(
                &self.hypergraphs[self.a],
                &self.hypergraphs[self.b],
                self.ring_operator.clone(),
                self.max_homology_dimension,
            );
        // due to the way distributed AM perform serialization and deserialization, we cannot just return the cospan directly
        // as it is generic over the ring element type, so we need to manually serialize it ourselfs, since we know at compile time the type
        let vec_bytes = bincode::serialize(&cospan).unwrap();
        let orig_cospan = bincode::deserialize::<Cospan<RingElement>>(&vec_bytes).expect("should be Cospan<RingElement>");
        vec_bytes
    }
}

#[lamellar::AmData]
//this is kind-of a ugly hack because we have to hardcode the ring element type as a bool
//do to how serializing and deserializing work with generics, hopefully we can eventually get rid of this
//(LocalAms do not have this constraint as they are not de/serialized)
struct CollectOnRootBoolRingElementAm { 
    root_data: LocalRwDarc<Vec<(usize,Vec<u8>)>>,
    // data: Vec<Cospan<bool>>,
    data: Vec<Vec<u8>>, //ugly ugly ugly
    pe: usize,
}

#[lamellar::am]
impl LamellarAM for CollectOnRootBoolRingElementAm {
    async fn exec(self) {
        let mut root_data = self.root_data.write().await;
        for (i, cospan) in self.data.iter().enumerate() {
            let offset = self.pe + lamellar::num_pes *i ;
            root_data.push((offset,cospan.clone()));
        }
    }
}

pub fn interval_decomposition_for_zigzag_of_hypgeraph_unions<RingOperator, RingElement>(
    hypergraphs: Vec<Vec<Vec<usize>>>,
    ring_operator: RingOperator,
    max_homology_dimension: usize,
    print_profiling_statistics:     bool,
) -> Vec<(usize, IntervalSubmoduleBasis<RingElement>)>
where
    RingOperator: Clone + Send + Sync + DivisionRing<RingElement> + std::fmt::Debug + 'static,
    RingElement: Clone + Send + Sync + Serialize + DeserializeOwned + std::fmt::Debug + Eq + Ord + 'static,
{
    let world = WORLD.clone();
    let num_pes = world.num_pes();
    let my_pe = world.my_pe();
    // compute number of vertices and arrows
    let n_hypergraphs = hypergraphs.len();

    // handle the edge case where the number of hypergraphs is 0 or 1
    // ========================================================================
    if n_hypergraphs == 0 {
        // if there are no hypergraphs then short circuit
        return Vec::with_capacity(0);
    } else if n_hypergraphs == 1 {
        let factored = factor_dowker_complex(
            hypergraphs[0].clone(),
            ring_operator.clone(),
            max_homology_dimension,
        );
        let mut submodules = Vec::with_capacity(factored.indices_harmonic().count());
        for basis_vector in factored.basis_harmonic_vec_of_vec().vec_of_vec() {
            let dimension = chain_dimension(basis_vector);
            let leftmost_vertex = 0;
            let submodule =
                IntervalSubmoduleBasis::new(leftmost_vertex, vec![basis_vector.clone()]);
            submodules.push((dimension, submodule));
        }
        return submodules;
    }

    // calculate the cospans
    // ========================================================================
    world.barrier();
    let start = Instant::now();

    let mut cospan_futures = vec![];
    let hypergraphs = Arc::new(hypergraphs);
    for (a, b) in (my_pe..hypergraphs.len())
        .tuple_windows()
        .step_by(num_pes)
    {
        // cospans.push(induced_cospan(
        //     hypergraph_a,
        //     hypergraph_b,
        //     ring_operator.clone(),
        //     max_homology_dimension,
        // ));

        // cospan_futures.push(world.spawn(async_induced_cospan(
        //     hypergraph_a.clone(),
        //     hypergraph_b.clone(),
        //     ring_operator.clone(),
        //     max_homology_dimension,
        // )));
        cospan_futures.push(world.exec_am_local(
            InducedCospanAm {
                hypergraphs: hypergraphs.clone(),
                a,
                b,
                ring_operator: ring_operator.clone(),
                max_homology_dimension,
                _phantom: std::marker::PhantomData,
            },
        ));
    }
    let root_data=LocalRwDarc::new(world.team(), vec![]).expect("Could not create root data"); // maybe we can overlap computation with setting up the darc
    let cospans = world.block_on_all(cospan_futures);
    world.exec_am_pe(0, CollectOnRootBoolRingElementAm {
        root_data: root_data.clone(),
        data: cospans,
        pe: my_pe
    }).block();
    world.barrier();
    let mut root_data_lock = root_data.blocking_write();

    root_data_lock
        .sort_by(
            |a, b| 
            {
                a.0.cmp(&b.0) 
            }
        );      
    //this is the dirty hack were we need to convert back from Vec<u8> to the Generic RingElement type
    let mut cospans = root_data_lock.drain(..).map(|x| bincode::deserialize::<Cospan<RingElement>>(&x.1).expect("should be Cospan<RingElement>")).collect::<Vec<_>>();

    world.barrier();
    let duration = start.elapsed(); // Stop the timer
    if print_profiling_statistics{ println!("Time elapsed to build cospans: {:?}", duration); }

    // return the decomposition
    // ========================================================================

    if world.my_pe() != 0 { cospans.clear() }
    cospan_pipeline::interval_decomposition_for_zigzag_of_cospans(cospans, ring_operator, print_profiling_statistics)
}

#[cfg(test)]
mod tests {

    #[test]
    fn doc_test_sideways_ladder() {
        use crate::algebra::rings::operator_structs::field_prime_order::BooleanFieldOperator;
        use crate::algebra::zigzag::hypergraph_pipeline::interval_decomposition_for_zigzag_of_hypgeraph_unions;
        use crate::topology::simplicial::from::relation::sideways_ladder_edges;

        // Construct length-2 sequence of hypergraphs of the following form. (in this
        // case there are actually just simple undirected graphs)
        //
        // hypergraph 0:
        //
        // 1 ---- 3 ---- 5
        // |      |      |
        // |      |      |
        // 0 ---- 2 ---- 4
        //
        // hypergraph 1:
        //
        //        3 ---- 5 ---- 7
        //        |      |      |
        //        |      |      |
        //        2 ---- 4 ---- 6
        //
        // --------------------------------------------------------------------

        let number_of_ladders = 2;
        let holes_per_ladder = 2;

        let mut hypergraphs = Vec::new();
        for offset_from_left in 0..number_of_ladders {
            hypergraphs.push(sideways_ladder_edges(offset_from_left, holes_per_ladder));
        }

        // Compute the barcode for the corresponding zigzag:
        //
        // hypergraph_0 ----> ( hypergraph_0  U  hypergraph_1 ) <---- hypergraph_1
        //
        // --------------------------------------------------------------------

        let ring_operator = BooleanFieldOperator::new();
        let max_homology_dimension = 2;
        let interval_decomposition = interval_decomposition_for_zigzag_of_hypgeraph_unions(
            hypergraphs,
            ring_operator,
            max_homology_dimension,
        );
        let barcode = interval_decomposition
            .iter()
            .map(|(dim, submodule)| (dim.clone(), submodule.interval_endpoints()))
            .collect::<Vec<_>>();

        // Verify that the barcode is correct
        // Note that each pair (a,b) represents a half-open interval of form [a,b)
        // The actual barcode is
        //
        //
        // dimension 0:
        //
        //      0 ---- 1 ---- 2
        //
        // dimension 1
        //
        //      0 ---- 1
        //      0 ---- 1 ---- 2
        //             1 ---- 2
        //
        // --------------------------------------------------------------------

        let barcode_ground_truth: Vec<(usize, (usize, usize))> =
            vec![(0, (0, 3)), (1, (0, 2)), (1, (0, 3)), (1, (1, 3))];

        assert_eq!(barcode, barcode_ground_truth);
    }

    #[test]
    fn test_sideways_ladder() {
        use crate::algebra::rings::operator_structs::field_prime_order::BooleanFieldOperator;
        use crate::algebra::zigzag::hypergraph_pipeline::interval_decomposition_for_zigzag_of_hypgeraph_unions;
        use crate::topology::simplicial::from::relation::sideways_ladder_edges;

        let number_of_ladders = 3;
        let holes_per_ladder = 2;

        let mut hypergraphs = Vec::new();
        for offset_from_left in 0..number_of_ladders {
            hypergraphs.push(sideways_ladder_edges(offset_from_left, holes_per_ladder));
        }

        // Compute the barcode for the corresponding zigzag:
        // hypergraph_0 ----> ( hypergraph_0  U  hypergraph_1 ) <---- hypergraph_1
        let ring_operator = BooleanFieldOperator::new();
        let max_homology_dimension = 2;
        let interval_decomposition = interval_decomposition_for_zigzag_of_hypgeraph_unions(
            hypergraphs,
            ring_operator,
            max_homology_dimension,
        );
        let barcode = interval_decomposition
            .iter()
            .map(|(dim, submodule)| (dim.clone(), submodule.interval_endpoints()))
            .collect::<Vec<_>>();

        // Verify that the barcode is correct
        // Note that each pair (a,b) represents a half-open interval of form [a,b)
        //
        // Here's a "diagram" of the ladders, where each * represents a hole
        //
        //  laddder 0        : * *
        //  laddder 1 (union): * * *
        //  laddder 2        :   * *
        //  laddder 3 (union):   * * *
        //  laddder 4        :     * *
        //
        // By in specting this we can see that the barcode should be
        //
        // Dimension 0:
        //      0 ---- 1 ---- 2 ---- 3 ---- 4
        //
        // Dimension 2:
        //      0 ---- 1
        //      0 ---- 1 ---- 2 ---- 3
        //             1 ---- 2 ---- 3 ---- 4
        //                           3 ---- 4

        let barcode_ground_truth = vec![
            (0, (0, 5)),
            (1, (0, 2)),
            (1, (0, 4)),
            (1, (1, 5)),
            (1, (3, 5)),
        ];

        assert_eq!(barcode, barcode_ground_truth);
    }

    /// Similar to preceding test, but for every hypergraph we add a hyperdge that contains all the odd vertices.
    #[test]
    fn test_sideways_ladder_with_top_blob() {
        use crate::algebra::rings::operator_structs::field_prime_order::BooleanFieldOperator;
        use crate::algebra::zigzag::hypergraph_pipeline::interval_decomposition_for_zigzag_of_hypgeraph_unions;
        use crate::topology::simplicial::from::relation::sideways_ladder_edges;
        use num::Integer;
        use std::collections::HashSet;

        let number_of_ladders = 3;
        let holes_per_ladder = 2;

        let mut hypergraphs = Vec::new();

        // build the ladder hypergraphs
        for offset_from_left in 0..number_of_ladders {
            hypergraphs.push(sideways_ladder_edges(offset_from_left, holes_per_ladder));
        }

        // add the blob
        let mut blob = HashSet::new();
        blob.extend(
            hypergraphs
                .iter()
                .map(|x| x.iter().map(|y| y.iter().cloned()).flatten())
                .flatten()
                .filter(|x| x.is_odd()),
        );
        let mut blob: Vec<_> = blob.into_iter().collect();
        blob.sort();

        for hypergraph in hypergraphs.iter_mut() {
            hypergraph.push(blob.clone());
        }

        // Compute the barcode for the corresponding zigzag:
        // hypergraph_0 ----> ( hypergraph_0  U  hypergraph_1 ) <---- hypergraph_1
        let ring_operator = BooleanFieldOperator::new();
        let max_homology_dimension = 2;
        let interval_decomposition = interval_decomposition_for_zigzag_of_hypgeraph_unions(
            hypergraphs,
            ring_operator,
            max_homology_dimension,
        );
        let barcode = interval_decomposition
            .iter()
            .map(|(dim, submodule)| (dim.clone(), submodule.interval_endpoints()))
            .collect::<Vec<_>>();

        // Verify that the barcode is correct
        // Note that each pair (a,b) represents a half-open interval of form [a,b)
        //
        // Here's a "diagram" of the ladders, where each * represents a hole
        //
        //  laddder 0        : * *
        //  laddder 1 (union): * * *
        //  laddder 2        :   * *
        //  laddder 3 (union):   * * *
        //  laddder 4        :     * *
        //
        // By in specting this we can see that the barcode should be
        //
        // Dimension 0:
        //      0 ---- 1 ---- 2 ---- 3 ---- 4
        //
        // Dimension 2:
        //      0 ---- 1
        //      0 ---- 1 ---- 2 ---- 3
        //             1 ---- 2 ---- 3 ---- 4
        //                           3 ---- 4

        let barcode_ground_truth = vec![
            (0, (0, 5)),
            (1, (0, 2)),
            (1, (0, 4)),
            (1, (1, 5)),
            (1, (3, 5)),
        ];

        assert_eq!(barcode, barcode_ground_truth);
    }
}
