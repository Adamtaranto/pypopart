# PyPopART Implementation TODO

## Project Overview
Pure Python implementation of PopART (Population Analysis with Reticulate Trees) for constructing and visualizing haplotype networks from DNA sequence data.

## Phase 1: Project Setup & Infrastructure ✅ COMPLETED

### 1.1 Development Environment ✅
- [x] Set up project structure with proper package layout
- [x] Create `pyproject.toml` with project metadata and dependencies
- [x] Configure development dependencies (pytest, ruff, mypy)
- [x] Set up pre-commit hooks for code quality
- [ ] Create virtual environment setup instructions

### 1.2 Testing Infrastructure ✅
- [x] Set up pytest configuration
- [x] Create test directory structure matching src layout
- [x] Set up test fixtures for sample sequence data
- [x] Configure code coverage reporting (pytest-cov)
- [ ] Create sample datasets for testing (small, medium, large)

### 1.3 CI/CD Pipeline
- [ ] Set up GitHub Actions for automated testing
- [ ] Configure multi-Python version testing (3.9, 3.10, 3.11, 3.12)
- [ ] Add linting checks to CI pipeline
- [ ] Set up automated coverage reporting
- [ ] Configure release automation

### 1.4 Documentation Framework
- [ ] Set up Sphinx for documentation generation
- [x] Create basic README.md with project description
- [ ] Set up documentation structure (installation, usage, API reference)
- [ ] Configure Read the Docs integration

## Phase 2: Core Data Structures ✅ COMPLETED

### 2.1 Sequence Representation ✅ COMPLETED (merged from copilot/vscode1763095267042)
- [x] Implement `Sequence` class for DNA sequences
  - [x] Store sequence ID, data, and metadata
  - [x] Support IUPAC nucleotide codes
  - [x] Handle ambiguous characters
  - [x] Implement sequence comparison methods
- [x] Add sequence validation methods
- [x] Implement sequence manipulation (reverse complement, translation)
- [x] Create unit tests for Sequence class (17 tests, 93% coverage)

### 2.2 Alignment Management ✅ COMPLETED (merged from copilot/vscode1763095267042)
- [x] Implement `Alignment` class for multiple sequence alignments
  - [x] Store collection of aligned sequences
  - [x] Validate alignment consistency (equal lengths)
  - [x] Access sequences by index or ID
  - [x] Calculate alignment statistics (length, gaps, conserved sites)
- [x] Implement alignment slicing and subsetting
- [x] Add methods for alignment quality checks
- [x] Create unit tests for Alignment class (25 tests, 97% coverage)

### 2.3 Haplotype Representation ✅ COMPLETED
- [x] Implement `Haplotype` class
  - [x] Store unique sequence variants
  - [x] Track frequency/count information
  - [x] Store associated sample IDs
  - [x] Handle population/group assignments
- [x] Implement haplotype identification from alignment
- [x] Add methods for haplotype diversity calculations
- [x] Create unit tests for Haplotype class (28 tests, 96% coverage)

### 2.4 Network Data Structure ✅ COMPLETED
- [x] Implement `HaplotypeNetwork` class using NetworkX
  - [x] Store nodes (haplotypes) with attributes
  - [x] Store edges (connections) with distances/weights
  - [x] Track median vectors (inferred ancestral nodes)
  - [x] Store network metadata and parameters
- [x] Add network manipulation methods (add/remove nodes/edges)
- [x] Implement network validation
- [x] Create unit tests for HaplotypeNetwork class (36 tests, 94% coverage)

## Phase 3: Distance Calculation & Metrics ✅ COMPLETED

### 3.1 Basic Distance Metrics ✅ COMPLETED
- [x] Implement Hamming distance calculator
- [x] Implement pairwise distance matrix calculation
- [x] Add support for handling gaps and ambiguous characters
- [x] Optimize distance calculations for large datasets
- [x] Create unit tests with known distance examples (28 tests, 91% coverage)

### 3.2 Evolutionary Distance Models ✅ COMPLETED
- [x] Implement Jukes-Cantor correction
- [x] Implement Kimura 2-parameter (K2P) distance
- [x] Implement Tamura-Nei distance
- [x] Add support for custom substitution models
- [x] Create unit tests for each model

### 3.3 Distance Matrix Management ✅ COMPLETED
- [x] Implement `DistanceMatrix` class
  - [x] Store pairwise distances efficiently
  - [x] Support both symmetric and asymmetric matrices
  - [x] Provide fast lookup methods
- [x] Add matrix export/import functionality (CSV format)
- [x] Implement matrix visualization (matplotlib heatmap)
- [x] Create unit tests for DistanceMatrix class (41 tests total, 90% coverage)

## Phase 4: Network Construction Algorithms ✅ COMPLETED

### 4.1 Minimum Spanning Tree (MST) ✅ COMPLETED
- [x] Implement MST algorithm using Prim's or Kruskal's
- [x] Add support for different distance metrics
- [x] Handle tied distances appropriately
- [x] Optimize for large datasets
- [x] Create comprehensive unit tests (13 tests, 95% coverage)
- [x] Add example usage documentation

### 4.2 Minimum Spanning Network (MSN) ✅ COMPLETED
- [x] Implement MSN algorithm
  - [x] Start with MST
  - [x] Add alternative connections at same distance
  - [x] Remove redundant edges
- [x] Handle network complexity management
- [x] Add visualization of alternative connections
- [x] Create unit tests with known examples (10 tests, 95% coverage)
- [x] Document algorithm parameters

### 4.3 TCS (Statistical Parsimony) ✅ COMPLETED
- [x] Implement parsimony probability calculation
- [x] Implement connection limit estimation (95% parsimony)
- [x] Build network using parsimony criterion
- [x] Handle ambiguous connections
- [x] Add support for missing data
- [x] Create unit tests with published examples (11 tests, 94% coverage)
- [x] Document statistical assumptions

### 4.4 Median-Joining Network ✅ COMPLETED
- [x] Implement median vector inference
  - [x] Generate potential median sequences
  - [x] Calculate median joining criterion
  - [x] Add median vectors to network
- [x] Implement network simplification
- [x] Add epsilon parameter for complexity control
- [x] Optimize computational efficiency
- [x] Create unit tests with known networks (14 tests, 59% coverage)
- [x] Document algorithm parameters and tuning

### 4.5 Algorithm Comparison Framework ✅ COMPLETED
- [x] Create unified interface for all algorithms (NetworkAlgorithm base class)
- [x] Implement algorithm selection and configuration
- [ ] Add benchmarking utilities
- [ ] Create comparative examples
- [ ] Document when to use each algorithm

## Phase 5: Network Analysis & Statistics ✅ COMPLETED

### 5.1 Frequency and Sample Data ✅ COMPLETED
- [x] Implement haplotype frequency calculations
- [x] Add population/group assignment tracking
- [x] Calculate per-population frequencies
- [x] Implement sample size normalization
- [ ] Create visualization of frequency distributions (deferred to Phase 6)

### 5.2 Network Statistics ✅ COMPLETED
- [x] Calculate network diameter
- [x] Compute clustering coefficients
- [x] Identify central haplotypes (degree, betweenness, closeness, eigenvector centrality)
- [x] Calculate reticulation index
- [x] Implement diversity metrics (nucleotide, haplotype, Shannon diversity)
- [x] Create summary statistics report
- [x] Create unit tests (29 tests, 91% coverage)

### 5.3 Population Genetics Measures ✅ COMPLETED
- [x] Implement Tajima's D
- [x] Calculate Fu's Fs
- [x] Compute pairwise FST
- [x] Calculate FST matrix for all population pairs
- [x] Implement AMOVA framework
- [x] Add mismatch distribution analysis
- [x] Create unit tests for each measure (21 tests, 97% coverage)

### 5.4 Network Topology Analysis ✅ COMPLETED
- [x] Identify star-like patterns (perfect and partial stars)
- [x] Detect network partitions (connected components)
- [x] Calculate node centrality measures (degree, betweenness, closeness, eigenvector)
- [x] Identify potential ancestral nodes (composite scoring)
- [x] Find hub nodes and bridges
- [x] Identify bottleneck nodes (articulation points)
- [x] Create topology summary reports (21 tests, 93% coverage)

## Phase 6: GUI and Visualization

### Core GUI

- [ ] Create Dash application structure
- [ ] Implement file upload component
- [ ] Create algorithm selection interface
- [ ] Build parameter controls
- [ ] Add compute button and progress indicator

### GUI Features

- [ ] Implement statistics panel
- [ ] Add alignment viewer
- [ ] Create export options
- [ ] Error handling and validation
- [ ] User testing and refinement

### 6.1 Static Network Plots
- [ ] Implement matplotlib-based network plotting
  - [ ] Node size proportional to frequency
  - [ ] Color nodes by population/group
  - [ ] Edge thickness by distance
  - [ ] Show median vectors distinctly
- [ ] Add customization options (colors, sizes, labels)
- [ ] Implement different layout algorithms
- [ ] Export high-resolution images (PNG, PDF, SVG)
- [ ] Create gallery of example plots

### 6.2 Interactive Visualization
- [ ] Implement Plotly-based interactive networks
  - [ ] Hover information for nodes and edges
  - [ ] Zoom and pan capabilities
  - [ ] Click to show sequence details
  - [ ] Toggle population visibility
- [ ] Export interactive HTML files
- [ ] Add filtering and search functionality
- [ ] Create interactive examples

### 6.4 Layout Algorithms
- [ ] Implement force-directed layout
- [ ] Add radial/circular layout option
- [ ] Create hierarchical layout
- [ ] Implement custom layout algorithms
- [ ] Allow manual node positioning
- [ ] Save and load layout configurations

### 6.5 Legend and Annotations
- [ ] Create informative legends
- [ ] Add scale bars for mutations
- [ ] Include summary statistics in plots
- [ ] Support custom annotations
- [ ] Create publication-ready figure templates

## Phase 7: File I/O & Data Import/Export

### 7.1 Sequence File Readers
- [ ] Implement FASTA file reader using Biopython
- [ ] Implement NEXUS file reader
- [ ] Implement PHYLIP file reader
- [ ] Add support for GenBank format
- [ ] Handle compressed files (gzip, zip)
- [ ] Create robust error handling and validation
- [ ] Add progress indicators for large files

### 7.2 Alignment Validation
- [ ] Check sequence lengths match
- [ ] Validate character sets
- [ ] Detect and report alignment issues
- [ ] Implement auto-correction options
- [ ] Create validation report

### 7.3 Network Export Formats
- [ ] Export to GraphML format
- [ ] Export to GML format
- [ ] Export to Cytoscape format
- [ ] Export to JSON for web applications
- [ ] Export to CSV/TSV for statistics
- [ ] Create format conversion utilities

### 7.4 PopART Compatibility
- [ ] Read PopART project files (.nex with traits)
- [ ] Export networks in PopART-compatible format
- [ ] Handle trait/metadata properly
- [ ] Create conversion utilities
- [ ] Test with existing PopART datasets

### 7.5 Metadata Management
- [ ] Implement trait/metadata file reading
- [ ] Support CSV metadata files
- [ ] Link metadata to sequences
- [ ] Validate metadata consistency
- [ ] Export metadata with networks

## Phase 8: User Interface & Documentation

### 8.1 Command-Line Interface
- [ ] Design CLI using Click or argparse
- [ ] Implement commands for each algorithm
- [ ] Add progress bars for long operations
- [ ] Create informative help messages
- [ ] Implement configuration file support
- [ ] Add verbose/quiet modes
- [ ] Create shell completion scripts

### 8.2 Python API
- [ ] Design clean, intuitive API
- [ ] Create high-level convenience functions
- [ ] Implement pipeline/workflow classes
- [ ] Add method chaining support
- [ ] Create comprehensive docstrings
- [ ] Write API reference documentation

### 8.3 Jupyter Notebook Integration
- [ ] Create example notebooks for each algorithm
- [ ] Implement rich display methods for networks
- [ ] Add interactive widgets for parameter tuning
- [ ] Create tutorial notebooks
- [ ] Set up notebook testing

### 8.4 Optional GUI (Future)
- [ ] Design GUI mockups
- [ ] Choose framework (PyQt6, Tkinter, or web-based)
- [ ] Implement file loading interface
- [ ] Create network visualization panel
- [ ] Add analysis parameter controls
- [ ] Implement result export functionality

### 8.5 Documentation
- [ ] Write comprehensive README
- [ ] Create installation guide
- [ ] Write user guide with examples
- [ ] Create algorithm comparison guide
- [ ] Document all API functions
- [ ] Create troubleshooting guide
- [ ] Write contribution guidelines
- [ ] Create FAQ section

### 8.6 Tutorials and Examples
- [ ] Create basic usage tutorial
- [ ] Write tutorials for each algorithm
- [ ] Create visualization customization guide
- [ ] Write population genetics analysis tutorial
- [ ] Create real-world case studies
- [ ] Add video tutorials (optional)

## Phase 9: Performance & Optimization

### 9.1 Profiling
- [ ] Profile code with realistic datasets
- [ ] Identify bottlenecks
- [ ] Create performance benchmarks
- [ ] Document performance characteristics

### 9.2 Optimization
- [ ] Optimize distance calculations with NumPy
- [ ] Implement parallel processing where appropriate
- [ ] Add caching for repeated calculations
- [ ] Optimize network algorithms
- [ ] Consider Numba for critical paths

### 9.3 Memory Management
- [ ] Profile memory usage
- [ ] Implement memory-efficient data structures
- [ ] Add streaming for large files
- [ ] Create memory usage documentation

## Phase 10: Validation & Testing

### 10.1 Algorithm Validation
- [ ] Test with published datasets
- [ ] Compare results with original PopART
- [ ] Validate against literature examples
- [ ] Create regression test suite
- [ ] Document any differences from PopART

### 10.2 Integration Testing
- [ ] Test complete workflows
- [ ] Test with diverse file formats
- [ ] Test with edge cases
- [ ] Test error handling
- [ ] Create end-to-end test suite

### 10.3 User Testing
- [ ] Conduct beta testing with researchers
- [ ] Gather feedback on usability
- [ ] Test on different platforms
- [ ] Create bug report template
- [ ] Implement feedback improvements

## Phase 11: Deployment & Distribution

### 11.1 Package Distribution
- [ ] Publish to PyPI
- [ ] Create conda package
- [ ] Set up automatic releases
- [ ] Create versioning strategy
- [ ] Write release notes template

### 11.2 Platform Support
- [ ] Test on Windows, macOS, Linux
- [ ] Create platform-specific installation guides
- [ ] Handle platform-specific issues
- [ ] Document platform requirements

### 11.3 Docker Support
- [ ] Create Dockerfile
- [ ] Publish to Docker Hub
- [ ] Create Docker usage documentation
- [ ] Add Docker Compose examples

## Dependencies to Consider

### Core Dependencies
- **NumPy**: Numerical operations and arrays
- **Biopython**: Sequence I/O and manipulation
- **NetworkX**: Graph/network data structures and algorithms
- **pandas**: Data manipulation and analysis

### Visualization
- **matplotlib**: Static plotting
- **plotly**: Interactive visualization
- **seaborn**: Statistical visualization

### Optional Dependencies
- **numba**: JIT compilation for performance
- **scipy**: Scientific computing and optimization
- **scikit-learn**: Machine learning algorithms (clustering, etc.)
- **geopandas**: Geographic data handling

### Development Tools
- **pytest**: Testing framework
- **black**: Code formatting
- **ruff**: Fast Python linter
- **mypy**: Type checking
- **sphinx**: Documentation generation

## Success Metrics

- [ ] All network algorithms produce correct results vs. known examples
- [ ] Performance handles datasets with 1000+ sequences
- [ ] Code coverage above 80%
- [ ] Comprehensive documentation with examples
- [ ] Successful validation against PopART outputs
- [ ] Positive user feedback from beta testers
- [ ] Published on PyPI with active downloads

## Notes

- Prioritize correctness over performance initially
- Maintain compatibility with PopART file formats where possible
- Design for extensibility (easy to add new algorithms)
- Follow Python best practices and PEP standards
- Consider eventual publication in bioinformatics journal
