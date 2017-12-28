import sys

"""
input:
C24882126       Cufflinks       gene    33      1196    .       +       .       ID=Gb_00001;
C24882126       Cufflinks       exon    33      116     .       +       0       Parent=Gb_00001;
C24882126       Cufflinks       intron  117     218     .       +       0       Parent=Gb_00001;
C24882126       Cufflinks       exon    219     460     .       +       0       Parent=Gb_00001;
C24882126       Cufflinks       intron  461     541     .       +       1       Parent=Gb_00001;
C24882126       Cufflinks       exon    542     863     .       +       1       Parent=Gb_00001;
C24882126       Cufflinks       intron  864     944     .       +       0       Parent=Gb_00001;
C24882126       Cufflinks       exon    945     1196    .       +       0       Parent=Gb_00001;

output:
C24882126       Cufflinks       gene    33      1196    .       +       .       ID=Gb_00001;
C24882126       Cufflinks       exon    219     460     .       +       0       ID=Gb_00001.exon.1;Parent=Gb_00001
C24882126       Cufflinks       intron  461     541     .       +       1       ID=Gb_00001.intron.1;Parent=Gb_00001
C24882126       Cufflinks       exon    542     863     .       +       1       ID=Gb_00001.exon.2;Parent=Gb_00001
C24882126       Cufflinks       intron  864     944     .       +       0       ID=Gb_00001.intron.2;Parent=Gb_00001
C24882126       Cufflinks       exon    945     1196    .       +       0       ID=Gb_00001.exon.3;Parent=Gb_00001
"""
# python input output

OUT = open(sys.argv[2], 'w')
exon_num = 0
intron_num = 0
tmp_lis1 = []
tmp_lis2 = []

with open(sys.argv[1], 'r') as f:
    for line in f:
        ls = line.split()
        feature = ls[8].split("=")[1].rstrip(";")
        if ls[2] == "gene":
            OUT.write("\t".join(ls) + "\n")
        elif ls[2] == "exon":
            if feature in tmp_lis1:
                exon_num += 1
                exon_id = feature + ".exon." + str(exon_num)
                OUT.write("\t".join(ls[:8]) + "\t" + "ID=" + exon_id + ";Parent=" + feature + "\n")
            else:
                tmp_lis1.append(feature)
                exon_num = 0
        elif ls[2] == "intron":
            if feature in tmp_lis2:
                intron_num += 1
                intron_id = feature + ".intron." + str(intron_num)
                OUT.write("\t".join(ls[:8]) + "\t" + "ID=" + intron_id + ";Parent=" + feature + "\n")
            else:
                tmp_lis2.append(feature)
                intron_num = 0

OUT.close()
