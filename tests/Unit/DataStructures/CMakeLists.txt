# Distributed under the MIT License.
# See LICENSE.txt for details.

set(DATA_STRUCTURES_TESTS
    DataStructures/Test_DataBox.cpp
    DataStructures/Test_DataVector.cpp
    DataStructures/Test_Index.cpp
    DataStructures/Test_IndexIterator.cpp
    DataStructures/Test_Matrix.cpp
    DataStructures/Test_SliceIterator.cpp
    DataStructures/Test_StripeIterator.cpp
    DataStructures/Test_Tensor.cpp
    DataStructures/Test_TensorExpressions.cpp
    DataStructures/Test_Variables.cpp
    )

set(SPECTRE_TESTS "${SPECTRE_TESTS};${DATA_STRUCTURES_TESTS}" PARENT_SCOPE)
