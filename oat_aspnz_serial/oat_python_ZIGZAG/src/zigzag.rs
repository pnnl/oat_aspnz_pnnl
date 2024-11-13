use std::time::Instant;

use oat_rust::algebra::zigzag::hypergraph_pipeline::interval_decomposition_for_zigzag_of_hypgeraph_unions_WITH_SPANS;
use oat_rust::algebra::{rings::operator_structs::field_prime_order::BooleanFieldOperator, zigzag::hypergraph_pipeline::interval_decomposition_for_zigzag_of_hypgeraph_unions};
use pyo3::types::IntoPyDict;
use pyo3::{pyfunction, types::PyDict, Py, PyAny, PyResult, Python};


use pyo3::prelude::*;






#[pyfunction]
pub fn barcode_for_union_zigzag_of_hypergraphs_z2_coefficients(
        hypergraphs:                Vec<Vec<Vec<usize>>>,
        max_homology_dimension:     usize,
    ) -> 
    PyResult<Vec< ( usize, (usize,usize) ) >>
{
    let ring_operator                                   =   BooleanFieldOperator::new();
    let interval_decomposition                          =   interval_decomposition_for_zigzag_of_hypgeraph_unions(hypergraphs, ring_operator, max_homology_dimension);
        let barcode                                     =   interval_decomposition.iter()
                                                                .map(| (dim, submodule)| ( dim.clone(), submodule.interval_endpoints() ) )
                                                                .collect::<Vec<_>>();
    Ok( barcode )    
}








/// Converts a sequence of hypergraphs into the corresponding zigzag homology barcodes, 
/// with cycle representatives.
///  
/// Homology is computed with coefficients in the two-element field.
/// 
/// The simplicial complexes in the construction are the downward-closures of the 
/// hypergraphs.
/// 
/// # Arguments
/// 
/// `hypergraphs`: a list of hypergraphs. Each hypergraph is formatted
/// as a list hyperedges, where each hyperedge is represented by a sorted list of integers. 
/// 
/// `max_homology_dimension`: the maximum dimension in which to compute homology. Homology
/// in dimensions 2 and above is often substantially more expensive in time and memory than
/// homology in dimension 1.
/// 
/// # Output
/// 
/// Given a sequence of hypergraphs `A, B, C, ..`, this method decomposes the zigzag
/// homology module for
/// 
/// A --> A u B <-- B --> B u C <-- ..
/// 
/// as a direct sum of interval modules. It then returns a Pandas DataFrame frame with
/// columns for
/// 
/// - `submodule uuid`: a unique integer assigned to each subdule, used as its "name"
/// - `dimension`: the homological dimension of the submodule
/// - `birth (inclusive)`: the first index at which the submodule is nonzero. For example,
///   if the feature first appears in hypergraph `A` this index will be zero, and if it 
///   first appears in hypergraph `A u B` then the index will be 0.5.
/// - `death (exclusive)`: the first index `i` such that the submodule remains zero at 
///   index `i` and all indices greater than `i`. For example, if the feature is present
///   in hypergraph `A` but not in hypergraph `A u B`, then the value of `death (exlusive)` 
///   will be 0.5
/// - `cycle representatives`: a Pandas DataFrame containing a list of cycle representatives
///   for an interval submodule. We think of the module as a representation of a quiver (aka 
///   directed graph) whose vertices are labeled 0, 0.5, 1, 1.5, .. as follows:
/// 
///   0 --> 0.5 <-- 1 --> 1.5 <-- .. 
///   
///   For each vertex in the quiver, the submodule contains a vector space of dimension 1
///   or zero. The `cycle representatives` DataFrame contains a basis vector (cycle 
///   representative) for each vertex where the subspace is nonzero.  The basis vector is
///   formatted as a list of simplices. This formatting represents a linear combination of 
///   simplices with coefficients in the two element field.
#[pyfunction]
pub fn homology_decomposition_for_union_zigzag_of_hypergraphs_z2_coefficients( 
        hypergraphs:                    Vec<Vec<Vec<usize>>>,
        max_homology_dimension:         usize,
        return_cycle_representatives:   bool,
        py:                             Python<'_>,        
    )  
    -> Py<PyAny> 
    {
        let start_omnia                                 =   Instant::now(); // Start the timer for formatting output                                                            

        let start_pandas                                 =   Instant::now();
        let pandas = py.import("pandas").ok().unwrap();          
        println!("Time to import pandas: {:?}", start_omnia.elapsed()); // Output the duration          


        let ring_operator                                   =   BooleanFieldOperator::new();
        let interval_decomposition                      =   interval_decomposition_for_zigzag_of_hypgeraph_unions(
                                                                hypergraphs, 
                                                                ring_operator, 
                                                                max_homology_dimension
                                                            );



        let start_format_output                         =   Instant::now(); // Start the timer for formatting output                                                            

        let mut homology_dimensions   =   Vec::with_capacity( interval_decomposition.len() );
        let mut cycle_representative_list_for_each_bar   =   Vec::with_capacity( interval_decomposition.len() );
        let mut birth   =   Vec::with_capacity( interval_decomposition.len() );
        let mut exclusive_death   = Vec::with_capacity( interval_decomposition.len() );


            // build the data columns        
        for ( dimension, submodule ) in interval_decomposition {


            let (b,d) = submodule.interval_endpoints();
            birth.push(b);
            exclusive_death.push(d);
            homology_dimensions.push( dimension );


            if return_cycle_representatives {
                let supporting_vertices: Vec<_>             =   ( b .. d ).collect();
                let cycle_representatives: Vec<_>           =   submodule.basis_vectors
                                                                    .into_iter()
                                                                    .map( 
                                                                        | chain |
                                                                        chain.into_iter().map( |(simplex,_)| simplex )
                                                                        .collect::<Vec<_>>()
                                                                    )
                                                                    .collect();
    
                let dict = PyDict::new(py);
                dict.set_item( "quiver vertex", supporting_vertices ).ok().unwrap();  
                dict.set_item( "cycle representative", cycle_representatives ).ok().unwrap();                
    
                let cycle_representative_list_for_single_bar: Py<PyAny> = pandas.call_method("DataFrame", ( dict, ), None)
                    .map(Into::into).ok().unwrap();   
    
                let kwarg = vec![("inplace", true)].into_py_dict(py);               
                cycle_representative_list_for_single_bar.call_method( py, "set_index", ( "quiver vertex", ), Some(kwarg)).ok().unwrap();                  
    
                cycle_representative_list_for_each_bar.push( cycle_representative_list_for_single_bar );  
            }                  

        }

        let birth           =   birth.into_iter().map(|x| x as f64 / 2.0 ).collect::<Vec<_>>();
        let exclusive_death           =   exclusive_death.into_iter().map(|x| x as f64 / 2.0 ).collect::<Vec<_>>();        

        let dict = PyDict::new(py);
        dict.set_item( "dimension", homology_dimensions ).ok().unwrap();          
        dict.set_item( "birth (inclusive)", birth ).ok().unwrap();  
        dict.set_item( "death (exclusive)", exclusive_death ).ok().unwrap();  

        if return_cycle_representatives{
            dict.set_item( "cycle representatives", cycle_representative_list_for_each_bar ).ok().unwrap();                
        }        
        

        

        let start_df                                 =   Instant::now();        
        let df: Py<PyAny>  = pandas.call_method("DataFrame", ( dict, ), None)
            .map(Into::into).ok().unwrap();
        println!("Time to create data frame: {:?}", start_df.elapsed()); // Output the duration          
        let index = df.getattr( py, "index", ).unwrap();
        index.setattr( py, "name", "submodule uuid", ).unwrap();
        println!("Time to create data frame + reindex: {:?}", start_df.elapsed()); // Output the duration                  

        println!("Time to format output: {:?}", start_format_output.elapsed()); // Output the duration
        println!("Time for everything but Python handoff: {:?}", start_omnia.elapsed()); // Output the duration                                                                    

        df
}    




















#[pyfunction]
pub fn homology_decomposition_for_union_zigzag_of_hypergraphs_z2_coefficients_WITH_SPANS( 
        hypergraphs:                    Vec<Vec<Vec<usize>>>,
        max_homology_dimension:         usize,
        py:                             Python<'_>,        
    )  
    -> Py<PyAny> 
    {
        let pandas = py.import("pandas").ok().unwrap();  


        let ring_operator                                   =   BooleanFieldOperator::new();
        let interval_decomposition                      =   interval_decomposition_for_zigzag_of_hypgeraph_unions_WITH_SPANS(
                                                                hypergraphs, 
                                                                ring_operator, 
                                                                max_homology_dimension
                                                            );


        let start_format_output                         =   Instant::now(); // Start the timer for formatting output
                                                    

        let mut homology_dimensions   =   Vec::with_capacity( interval_decomposition.len() );
        let mut cycle_representative_list_for_each_bar   =   Vec::with_capacity( interval_decomposition.len() );
        let mut birth   =   Vec::with_capacity( interval_decomposition.len() );
        let mut exclusive_death   = Vec::with_capacity( interval_decomposition.len() );


            // build the data columns        
        for ( dimension, submodule ) in interval_decomposition {


            let (b,d) = submodule.interval_endpoints();
            birth.push(b);
            exclusive_death.push(d);
            homology_dimensions.push( dimension );


            let supporting_vertices: Vec<_>             =   ( b .. d ).collect();
            let cycle_representatives: Vec<_>           =   submodule.basis_vectors
                                                                .into_iter()
                                                                .map( 
                                                                    | chain |
                                                                    chain.into_iter().map( |(simplex,_)| simplex )
                                                                    .collect::<Vec<_>>()
                                                                )
                                                                .collect();

            let dict = PyDict::new(py);
            dict.set_item( "quiver vertex", supporting_vertices ).ok().unwrap();  
            dict.set_item( "cycle representative", cycle_representatives ).ok().unwrap();                

            let cycle_representative_list_for_single_bar: Py<PyAny> = pandas.call_method("DataFrame", ( dict, ), None)
                .map(Into::into).ok().unwrap();   

            let kwarg = vec![("inplace", true)].into_py_dict(py);               
            cycle_representative_list_for_single_bar.call_method( py, "set_index", ( "quiver vertex", ), Some(kwarg)).ok().unwrap();                  

            cycle_representative_list_for_each_bar.push( cycle_representative_list_for_single_bar );                    

        }

        let birth           =   birth.into_iter().map(|x| x as f64 / 2.0 ).collect::<Vec<_>>();
        let exclusive_death           =   exclusive_death.into_iter().map(|x| x as f64 / 2.0 ).collect::<Vec<_>>();        

        let dict = PyDict::new(py);
        dict.set_item( "dimension", homology_dimensions ).ok().unwrap();          
        dict.set_item( "birth (inclusive)", birth ).ok().unwrap();  
        dict.set_item( "death (exclusive)", exclusive_death ).ok().unwrap();          
        dict.set_item( "cycle representatives", cycle_representative_list_for_each_bar ).ok().unwrap();                
        
        let df: Py<PyAny>  = pandas.call_method("DataFrame", ( dict, ), None)
            .map(Into::into).ok().unwrap();
        let index = df.getattr( py, "index", ).unwrap();
        index.setattr( py, "name", "submodule uuid", ).unwrap();

        println!("Time to format output: {:?}", start_format_output.elapsed()); // Output the duration                                                            

        df
}    
























#[pyfunction]
pub fn homology_decomposition_for_union_zigzag_of_hypergraphs_z2_coefficients_WITH_SPANS_barcode_only( 
        hypergraphs:                    Vec<Vec<Vec<usize>>>,
        max_homology_dimension:         usize,
        py:                             Python<'_>,        
    )  
    -> Py<PyAny> 
    {
        let pandas = py.import("pandas").ok().unwrap();  


        let ring_operator                                   =   BooleanFieldOperator::new();
        let interval_decomposition                      =   interval_decomposition_for_zigzag_of_hypgeraph_unions_WITH_SPANS(
                                                                hypergraphs, 
                                                                ring_operator, 
                                                                max_homology_dimension
                                                            );


        let start_format_output                         =   Instant::now(); // Start the timer for formatting output
                                                    

        let mut homology_dimensions   =   Vec::with_capacity( interval_decomposition.len() );
        let mut birth   =   Vec::with_capacity( interval_decomposition.len() );
        let mut exclusive_death   = Vec::with_capacity( interval_decomposition.len() );


            // build the data columns        
        for ( dimension, submodule ) in interval_decomposition {


            let (b,d) = submodule.interval_endpoints();
            birth.push(b);
            exclusive_death.push(d);
            homology_dimensions.push( dimension );
            
        }

        let birth           =   birth.into_iter().map(|x| x as f64 / 2.0 ).collect::<Vec<_>>();
        let exclusive_death           =   exclusive_death.into_iter().map(|x| x as f64 / 2.0 ).collect::<Vec<_>>();        

        let dict = PyDict::new(py);
        dict.set_item( "dimension", homology_dimensions ).ok().unwrap();          
        dict.set_item( "birth (inclusive)", birth ).ok().unwrap();  
        dict.set_item( "death (exclusive)", exclusive_death ).ok().unwrap();                        
        
        let df: Py<PyAny>  = pandas.call_method("DataFrame", ( dict, ), None)
            .map(Into::into).ok().unwrap();
        let index = df.getattr( py, "index", ).unwrap();
        index.setattr( py, "name", "submodule uuid", ).unwrap();

        println!("Time to format output: {:?}", start_format_output.elapsed()); // Output the duration                                                            

        df
}    