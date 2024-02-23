## RTRex

This repositary contains the codes for running the RTRex (Regularly Triangle Rich Extractor) for undirected graphs.

## Input Data

The graph file must have the following structure:

* An initial line containing: “<number of vertices> <number of edges>”
* For each edge in the graph will contain an additional line of the type: “<src vertex id> <dest vertex it>”

The graph is assumed to be undirected. The vertices Id’s must go from 0 to <number of vertices>-1. 
Please refer to the Spectral-Triadic/Example/small-test.txt as an example

## Compile the Code

In RTex/ run:
`make`

## Run the RTRex code

In ‘src/RTRex/clustering/’ run:
`./RTRex <graph_dir> <name> <epsilon> <mode> <additional parameters>`

Parameters:

* <graph_dir>: Address of the file containing the graph (see Input Requirements)
* <name>: Name for the output files
* <epsilon>: Value of epsilon for the algorithm. Recommended between 0.1-0.3 (I used 0.1 mostly)
* <mode>: Different modes, we will be adding more modes here for different variants:
    * ‘s’: Simple or Standard mode, it generates non-overlapping decompositions. (Default if not specified)
    * ‘e’: Expanded, it runs the same decomposition than the Standard mode but it will also run the code to add the Singleton vertices to some clusters. It requires an additional parameter
* <threshold>: that specifies how many neighbors are required for the isolated vertices. I used 20 for eg_12_200. Best will depend on the average degree of the graph.

Examples:
`./RTRex ../Example/small-test.txt example 0.1 s`
`./RTex ../Example/small-test.txt example 0.1 e 20`

## Output

The code generates a single output file in the specified address. The file contains 2 lines for each cluster:

* First line is of the type: “Cluster id: 0, C: 2, n_vertices: 4, n_edges: 6, parent: 1”
* Second line contains a list of all the vertices in the correspondent cluster.

## Security

See [CONTRIBUTING](CONTRIBUTING.md#security-issue-notifications) for more information.

## License

This project is licensed under the Apache-2.0 License.

