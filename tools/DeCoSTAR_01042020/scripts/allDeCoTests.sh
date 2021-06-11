echo "********************* base tests *********************"
echo "--------------------- DeCo without transfer ---------------------"
bin/DeCoSTAR parameter.file=tests/testDeCo.param.txt
echo "--------------------- DeCo with transfer ---------------------"
bin/DeCoSTAR parameter.file=tests/testDeCoLT.param.txt

echo "********************* various options and algorithms *********************"
echo "--------------------- DeCo with transfer but without time slices ---------------------"
bin/DeCoSTAR parameter.file=tests/testDeCo_LT_undated.param.txt
echo "--------------------- DeCo with boltzmann sampling without transfer ---------------------"
bin/DeCoSTAR parameter.file=tests/testDeCo_boltzmann.param.txt
echo "--------------------- DeCo with boltzmann sampling with transfer ---------------------"
bin/DeCoSTAR parameter.file=tests/testDeCoLT_boltzmann.param.txt
echo "--------------------- DeCo with transfer and bounded time slices ---------------------"
bin/DeCoSTAR parameter.file=tests/testDeCoLT_BTS.param.txt
echo "--------------------- DeCo with transfer ---------------------"
bin/DeCoSTAR parameter.file=tests/testDeCo_LT_dated.param.txt
echo "--------------------- DeCo with transfer and fancy parameters ---------------------"
bin/DeCoSTAR parameter.file=tests/testDeCoLTfancy.param.txt
echo "--------------------- DeCo with different adjacency transmission rules ------------"
bin/DeCoSTAR parameter.files=tests/testDeCo.TransmissionAll.param.txt

echo "********************* reconciled trees as input *********************"
echo "--------------------- DeCo with transfer and reconciled trees ---------------------"
bin/DeCoSTAR parameter.file=tests/testDeCoLT_recInput.param.txt
echo "--------------------- DeCo with transfer, bounded time slices and reconciled trees ---------------------"
bin/DeCoSTAR parameter.file=tests/testDeCoLT_recInput_BTS.param.txt
echo "--------------------- DeCo without transfer and reconciled trees ---------------------"
bin/DeCoSTAR parameter.file=tests/testDeCoLT_recInput_NOTR.param.txt
echo "--------------------- DeCo with transfer, without time slices and reconciled trees ---------------------"
bin/DeCoSTAR parameter.file=tests/testDeCoLT_recInput_undated.param.txt

echo "********************* ARt-DeCo : discovering adjacencies in badly assembled extant genomes *********************"
bin/DeCoSTAR parameter.file=tests/Anopheles_dataset/testDeCo_ARt-DeCo.param.txt

echo "********************* ADseq : confidence score on scaffolding derived adjacencies *********************"
echo "--------------------- basic version ---------------------"
bin/DeCoSTAR parameter.file=tests/testADseq.param.txt
echo "--------------------- version with boltzmann sampling ---------------------"
bin/DeCoSTAR parameter.file=tests/testADseq_boltzmann.param.txt


echo "********************* ale input - longer test (a few seconds...) *********************"
bin/DeCoSTAR parameter.file=tests/testDeCoLT_aleinput.param.txt

echo "********************* Loss aware : free adjacency gains between the neighbours of a loss under certain conditions ************************"
./bin/DeCoSTAR parameter.file=tests/LossAware.test.param.txt


echo "tests finished. You may see the results in tests/testResults/ in various folders. You can compare with the content of expected/ ."
echo "expect some difference between adjacencies histories because of the random nature of the process (co-optimal solutions, or boltzmann sampling)."
