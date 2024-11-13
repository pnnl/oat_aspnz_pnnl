use oat_rust::algebra::{
    rings::operator_structs::field_prime_order::BooleanFieldOperator,
    zigzag::hypergraph_pipeline::interval_decomposition_for_zigzag_of_hypgeraph_unions,
};
use pyo3::types::IntoPyDict;
use pyo3::{pyfunction, types::PyDict, Py, PyAny, PyResult, Python};

use pyo3::prelude::*;

/// Converts a sequence of hypergraphs into its corresponding 0 and 1 dimensional barcodes.
///
/// Concretely, we take a sequence of hypergraphs `A, B, C, ..` and decompose the zigzag module for
/// `A --> A u B <-- B --> B u C <-- ..`.
///
/// The coefficient field is the 2-element field.
///
/// The simplicial complexes in the construction are the downward-closures fo the hypergraphs.
///
/// # Input requirements
///
/// The input to this function is a list of list-of-lists-of-integers. Each list-of-list-of-integers
/// represents a hypergraph. Each list of integers it contains must be sorted.
#[pyfunction]
pub fn barcode_for_union_zigzag_of_hypergraphs_z2_coefficients(
    hypergraphs: Vec<Vec<Vec<usize>>>,
    max_homology_dimension: usize,
    print_profiling_statistics: bool,
) -> PyResult<Vec<(usize, (usize, usize))>> {
    let ring_operator = BooleanFieldOperator::new();
    let interval_decomposition = interval_decomposition_for_zigzag_of_hypgeraph_unions(
        hypergraphs,
        ring_operator,
        max_homology_dimension,
        print_profiling_statistics,
    );
    let barcode = interval_decomposition
        .iter()
        .map(|(dim, submodule)| (dim.clone(), submodule.interval_endpoints()))
        .collect::<Vec<_>>();
    Ok(barcode)
}

// /// A function that creates and returns a Pandas DataFrame
// #[pyfunction]
// fn create_dataframe(py: Python) -> PyResult<PyObject> {
//     // Import pandas module
//     let pandas = py.import("pandas")?;

//     // Define some data as a dictionary of lists
//     let data = PyDict::new(py);
//     data.set_item("column1", vec![1, 2, 3])?;
//     data.set_item("column2", vec![4, 5, 6])?;

//     // Create the DataFrame using pandas.DataFrame
//     let df = pandas.get("DataFrame")?.call1((data,))?;

//     // Return the DataFrame as a PyObject
//     Ok(df.into())
// }

// /// A Python module implemented in Rust.
// #[pymodule]
// fn my_rust_module(py: Python, m: &PyModule) -> PyResult<()> {
//     m.add_function(wrap_pyfunction!(create_dataframe, m)?)?;
//     Ok(())
// }

#[pyfunction]
pub fn homology_decomposition_for_union_zigzag_of_hypergraphs_z2_coefficients(
    hypergraphs: Vec<Vec<Vec<usize>>>,
    max_homology_dimension: usize,
    return_cycle_representatives: bool,
    print_profiling_statistics: bool,
    py: Python<'_>,
) -> Py<PyAny> {
    let pandas = py.import("pandas").ok().unwrap();

    let ring_operator = BooleanFieldOperator::new();
    let interval_decomposition = interval_decomposition_for_zigzag_of_hypgeraph_unions(
        hypergraphs,
        ring_operator,
        max_homology_dimension,
        print_profiling_statistics,
    );

    let mut homology_dimensions = Vec::with_capacity(interval_decomposition.len());
    let mut cycle_representative_list_for_each_bar =
        Vec::with_capacity(interval_decomposition.len());
    let mut birth = Vec::with_capacity(interval_decomposition.len());
    let mut exclusive_death = Vec::with_capacity(interval_decomposition.len());

    // build the data columns
    for (dimension, submodule) in interval_decomposition {
        let (b, d) = submodule.interval_endpoints();
        birth.push(b);
        exclusive_death.push(d);
        homology_dimensions.push(dimension);

        let supporting_vertices: Vec<_> = (b..d).collect();

        if return_cycle_representatives {
            let cycle_representatives: Vec<_> = submodule
            .basis_vectors
            .into_iter()
            .map(|chain| {
                chain
                    .into_iter()
                    .map(|(simplex, _)| simplex)
                    .collect::<Vec<_>>()
            })
            .collect();

            let dict = PyDict::new(py);
            dict.set_item("quiver vertex", supporting_vertices)
                .ok()
                .unwrap();
            dict.set_item("cycle representative", cycle_representatives)
                .ok()
                .unwrap();

            let cycle_representative_list_for_single_bar: Py<PyAny> = pandas
                .call_method("DataFrame", (dict,), None)
                .map(Into::into)
                .ok()
                .unwrap();

            let kwarg = vec![("inplace", true)].into_py_dict(py);
            cycle_representative_list_for_single_bar
                .call_method(py, "set_index", ("quiver vertex",), Some(kwarg))
                .ok()
                .unwrap();

            cycle_representative_list_for_each_bar.push(cycle_representative_list_for_single_bar);
        }


        // // basis vectors
        // let enumerated_basis_vectors                               =   decomposition.basis_vectors_for_submodule(&single_bar_ledger);

        // let mut vertices    =   Vec::new();
        // let mut basis_vectors: Vec<Vec<usize>>   =   Vec::new();

        // for (vertex,chain) in enumerated_basis_vectors {
        //     vertices.push(vertex);
        //     let simplices = chain.into_iter().map(|(x,y)| x ).collect::<Vec<_>>();
        //     basis_vectors.push( simplices )
        // }

        // let dict = PyDict::new(py);
        // dict.set_item( "vertex", vertices ).ok().unwrap();
        // dict.set_item( "cycle representative", basis_vectors ).ok().unwrap();
        // let cycle_representative_list_for_single_bar: Py<PyAny> = pandas.call_method("DataFrame", ( dict, ), None)
        //     .map(Into::into).ok().unwrap();
        // cycle_representative_list_for_each_bar.push( cycle_representative_list_for_single_bar );
    }

    let dict = PyDict::new(py);
    dict.set_item("dimension", homology_dimensions)
        .ok()
        .unwrap();
    dict.set_item("birth (inclusive)", birth).ok().unwrap();
    dict.set_item("death (exclusive)", exclusive_death)
        .ok()
        .unwrap();

    if return_cycle_representatives {

        dict.set_item(
            "cycle representatives",
            cycle_representative_list_for_each_bar,
        )
        .ok()
        .unwrap();
    }    


    let df: Py<PyAny> = pandas
        .call_method("DataFrame", (dict,), None)
        .map(Into::into)
        .ok()
        .unwrap();
    let index = df.getattr(py, "index").unwrap();
    index.setattr(py, "name", "submodule uuid").unwrap();
    df
}
