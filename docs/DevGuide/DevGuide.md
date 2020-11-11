\cond NEVER
Distributed under the MIT License.
See LICENSE.txt for details.
\endcond
# Developer's Guide {#dev_guide}

### Developing and Improving Executables
- \ref spectre_build_system "Build system" and how to add dependencies,
  unit tests, and executables.
- \ref dev_guide_creating_executables "Executables and how to add them"
- \ref tutorials_parallel - A series of tutorials demonstrating how to
  write a parallel executable, and explaining some of the
  metaprogramming that turns user-provided code into a SpECTRE
  parallel executable
- \ref dev_guide_option_parsing "Option parsing" to get options from input files
- \ref dev_guide_importing "Importing" data from files
- \ref profiling_with_projections "Profiling With Charm++ Projections" and PAPI
  for optimizing performance
- \ref spectre_writing_python_bindings "Writing Python Bindings" to use
  SpECTRE C++ classes and functions from within python.
- \ref implementing_vectors "Implementing SpECTRE vectors" a quick how-to for
  making new generalizations of DataVectors
- \ref compiler_and_linker_errors "How to parse linker and compiler errors"

### Having your Contributions Merged into SpECTRE
- \ref writing_good_dox "Writing good documentation" is key for long term
  maintainability of the project.
- \ref writing_unit_tests "Writing Unit Tests" to catch bugs and to make
  sure future changes don't cause your code to regress.
- \ref travis_guide "Travis CI" is used to test every pull request.
- \ref code_review_guide "Code review guidelines." All code merged into
  develop must follow these requirements.

### General SpECTRE Terminology
Terms with SpECTRE-specific meanings are defined here.
- \ref domain_concepts "Domain Concepts" used throughout the code are defined
  here for reference.

### Template Metaprogramming (TMP)
Explanations for TMP concepts and patterns known to the greater C++ community
can be found here.
- \ref sfinae "SFINAE"
- \ref brigand "The Brigand TMP library"

### Foundational Concepts in SpECTRE
Designed to give the reader an introduction to SpECTRE's most recurring
concepts and patterns.
- \ref databox_foundations "Towards SpECTRE's DataBox"
- \ref protocols "Protocols: metaprogramming interfaces"

### Technical Documentation for Fluent Developers
Assumes a thorough familiarity and fluency in SpECTRE's usage of TMP.
- \ref DataBoxGroup "DataBox"
- \ref observers_infrastructure_dev_guide "Observers infrastructure"
- \ref dev_guide_parallelization_foundations "Parallelization infrastructure"
  components and usage

### CoordinateMap Guide
Methods for creating custom coordinate maps are discussed here.
- \ref redistributing_gridpoints "Methods for redistributing gridpoints"
