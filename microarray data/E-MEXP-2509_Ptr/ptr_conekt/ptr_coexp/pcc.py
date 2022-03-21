import numpy as np
import sys

filename = 'ptr_modified.txt'
output = 'ptr.PCC'
mcl_output = 'ptr.forMCL'

# Read Matrix and store nominators and denominators
with open(filename, 'r') as fin:
    nominators, denominators, genes = [], [], []

    header = fin.readline()
    size = len(header.strip().split('\t'))

    for line in fin:
            parts = line.rstrip().split("\t")

            if size != len(parts):
                print("Warning! Unequal number of columns found in line:\n%s.\nExpression matrix corrupt. Aborting!\n" % line, file=sys.stderr)
                quit()

            if len(parts) == size:
                temp = []
                for j in range(1, len(parts)):
                    try:
                        temp.append(float(parts[j]))
                    except ValueError:
                        print("Warning! Non-number character found in line:\n%s.\nExpression matrix corrupt. Aborting!\n" % line, file=sys.stderr)
                        quit()

                row_values = np.array(temp)
                nomi = row_values-(sum(row_values)/len(row_values))
                denomi = np.sqrt(sum(nomi**2))

                if denomi != 0.0:
                    nominators.append(nomi)
                    denominators.append(denomi)
                    genes.append(parts[0])

nominators = np.array(nominators)
denominators = np.array(denominators)

# Calculate PCC and write output
with open(output, 'w') as fout, open(mcl_output, 'w') as mcl_out:
    print("Database OK.\nCalculating Pearson Correlation Coefficient and ranks.\n")
    for i, (nom, denom, gene) in enumerate(zip(nominators, denominators, genes), start=1):
        print("Calculated PCC values for sequence:%s, %d out of %d." % (gene, i, len(nominators)))

        nominator = np.dot(nominators, nom)
        denominator = np.dot(denominators, denom)
        pcc_values = nominator/denominator

        data = [{'score': p,
                 'gene': g,
                 'string': g + '(' + str(p) + ')'} for g, p in zip(genes, pcc_values) if g != gene]

        # sort by absolute pcc value
        data.sort(key=lambda x: x['score'], reverse=True)

        # get top 1000 genes and write them
        subset = data[:1000]

        fout.writelines(gene + ": " + '\t'.join([s['string'] for s in subset]) + "\n")

        for s in subset:
            # Keep scores > 0.7 substract 0.7 from result to remap values to [0,0.3] as this is important for mcl
            if s['score'] > 0.7:
                print(gene, s['gene'], s['score'] - 0.7, sep='\t', file=mcl_out)

print("PCCs calculated and saved as %s and %s." % (output, mcl_output), file=sys.stderr)