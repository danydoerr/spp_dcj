INSTALL_DIR=$1

mkdir -p ${INSTALL_DIR}

# Bio++ core
cd bpp-core/
cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR}
make install
cd ../

# Bio++ seq
cd bpp-seq/
cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR}
make install
cd ../

# Bio++ phyl
cd bpp-phyl/
cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR}
make install
cd ../
