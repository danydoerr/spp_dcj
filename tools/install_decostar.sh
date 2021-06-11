BPP_DIR=$1

cd DeCoSTAR_01042020/
cp makefile makefile_v0
sed -i 's|BPP_INCLUDE=/include|BPP_INCLUDE='${BPP_DIR}'/include|g' makefile
sed -i 's|BPP_LIB=/lib|BPP_LIB='${BPP_DIR}'/lib|g' makefile
# We assume the path to boost libraries do not need to be configured
make bin/DeCoSTAR

