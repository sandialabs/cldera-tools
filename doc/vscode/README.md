# VSCode Setup Instructions
 - Install vscode: <https://code.visualstudio.com/download>
 - For remote vscode, follow these instructions: <https://code.visualstudio.com/docs/remote/ssh>

## E3SM Setup
 - Install the "Modern Fortran" and "FORTRAN IntelliSense" vscode extensions
 - Install "fortls". For a local install, use: `pip3 install fortran-language-server --user`
 - Copy `.fortls` and `e3sm-workspace.code-workspace` into the E3SM main directory
 - Edit `e3sm-workspace.code-workspace` to use your own paths and preferences
 - Edit `.fortls` to point to the directories you would like to parse through. More instructions for editting this file can be found here: <https://github.com/hansec/fortran-language-server>

## CLDERA-tools Setup
 - Install the "C/C++" and "clangd" vscode extensions for C++ parsing
 - Install "clangd" or on mappy use `/home/jwatkin/bin/clangd`
 - Create a `compile_commands.json` file by building cldera-tools with the following cmake configuration:
```
MPI_CXX_INCLUDE_DIRS=$(mpicxx -show | cut -d ' ' -f 2)
MPI_CC_INCLUDE_DIRS=$(mpicc -show | cut -d ' ' -f 2)
cmake \
...
  -D CMAKE_C_FLAGS="${MPI_CC_INCLUDE_DIRS}" \
  -D CMAKE_CXX_FLAGS="{$MPI_CXX_INCLUDE_DIRS}" \
  -D CMAKE_EXPORT_COMPILE_COMMANDS=ON \
...
```
 - Install the "Modern Fortran" and "FORTRAN IntelliSense" vscode extensions for FORTRAN parsing
 - Install "fortls". For a local install, use: `pip3 install fortran-language-server --user`
 - Copy `cldera-tools-workspace.code-workspace` into the cldera-tools main directory
 - Edit `cldera-tools-workspace.code-workspace` to use your own paths and preferences

