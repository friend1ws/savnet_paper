#! /usr/bin/env python

from genomon_splicing_mutation.utils import check_significance

import numpy

N = 20
lambda1 = 3
lambda2 = 0.3
output_dir = "../output/"

hout = open(output_dir + "sim1.input.txt", 'w') 
for k in range(1, 6):
    # active
    for n in range(1000):

        sp_count = []
        for j in range(k):
            active = numpy.random.poisson(lambda1, 1)
            inactive = numpy.random.poisson(lambda2, N - 1)
            sp_count.append(','.join([str(x) for x in numpy.r_[active, inactive]]))

        test_name = "active_" + str(k) + "_" + str(n)
        mutation_status = "1:1"
        link_status = ';'.join(["1," + str(x) for x in range(1, k + 1)])

        print >> hout, test_name + '\t' + mutation_status + '\t' + ';'.join(sp_count) + '\t' + link_status

    # inactive

    for n in range(1000):

        sp_count = []
        for j in range(k):
            active = numpy.random.poisson(lambda2, 1)
            inactive = numpy.random.poisson(lambda2, N - 1)
            sp_count.append(','.join([str(x) for x in numpy.r_[active, inactive]]))

        test_name = "inactive_" + str(k) + "_" + str(n)
        mutation_status = "1:1"
        link_status = ';'.join(["1," + str(x) for x in range(1, k + 1)])

        print >> hout, test_name + '\t' + mutation_status + '\t' + ';'.join(sp_count) + '\t' + link_status

hout.close()

weight_vector = [1.0] * 20
alpha0 = 1.0
beta0 = 1.0
alpha1 = 1.0
beta1 = 0.01

check_significance(output_dir + "sim1.input.txt", output_dir + "sim1.output.tmp.txt", weight_vector, 5.0, alpha0, beta0, alpha1, beta1)


hout = open(output_dir + "sim1.output.txt", 'w')
print >> hout, "Is_Active" + '\t' + "Splicing_Num" + '\t' + "Trial" + '\t' + "BF"
with open(output_dir + "sim1.output.tmp.txt", 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        FF = F[0].split('_')
        print >> hout, FF[0] + '\t' + FF[1] + '\t' + FF[2] + '\t' + F[5]

hout.close()


