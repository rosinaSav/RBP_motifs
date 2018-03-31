'''
Author: Rosina Savisaar.
Compile a set of motifs putatively recognized by RNA-binding proteins.
'''

from housekeeping import flatten, list_to_dict, overlap, parse_arguments, print_elements
import matplotlib.pyplot as plt
import nucleotide_comp as nc
import numpy as np
import os
import plotting
import re
import read_and_write as rw

def main():
    '''
    Read in a series of input files on the sequence specificities of RBPs,
    filter the data and write a set of motifs for each RBP.
    Arguments (see Methods for further details on the input data files):
    upper_threshold, lower_threshold: the longest and shortest a motif is allowed to be, respectively
    RBPDB_experiments: path to RBPDB experiments file
    RBPDB proteins: path to RBPDB proteins file
    RBPDB_PWMs: path to file containing RBPDB PWM identifier to RBP mapping
    pwm_dir: path to directory containing RBPDB PWMs
    RBPmap_PSSMs: path to directory containing RBPmap PSSMs
    SFmap_proteins: path to file containing motifs from SFmap
    RNAcompete_information: path to summary file from CIS-BP RNA
    RNAcompete_PWMs: path to directory containing CIS-BP RNA PWMs
    final_motifs_file_name: name for output file
    plot_name: file for plot displaying the distribution of motif set sizes
    species: the species for which motifs are required
    '''

    description = "Compile a set of motifs putatively recognized by RNA-binding proteins."
    args = parse_arguments(description, ["upper_threshold", "lower_threshold", "RBPDB_experiments", "RBPDB_proteins", "RBPDB_PWMs", "pwm_dir", "RBPmap_PSSMs", "SFmap_proteins", "RNAcompete_information", "RNAcompete_PWMs", "final_motifs_file_name", "plot_name", "species"], ints = [0, 1])
    [upper_threshold, lower_threshold, RBPDB_experiments, RBPDB_proteins, RBPDB_PWMs, pwm_dir, RBPmap_PSSMs, SFmap_proteins, RNAcompete_information, RNAcompete_PWMs, final_motifs_file_name, plot_name, species] = [args.upper_threshold, args.lower_threshold, args.RBPDB_experiments, args.RBPDB_proteins, args.RBPDB_PWMs, args.pwm_dir, args.RBPmap_PSSMs, args.SFmap_proteins, args.RNAcompete_information, args.RNAcompete_PWMs, args.final_motifs_file_name, args.plot_name, args.species]

    db_fields = rw.read_many_fields(RBPDB_experiments, ",")
    db_fields = db_fields[1:]
    print("There are {0} RBPDB experiments.".format(len(db_fields)))
    db_proteins = rw.read_many_fields(RBPDB_proteins, ",")
    #species is "Homo sapiens" or "Mus musculus"
    db_proteins = [i for i in db_proteins if i[6] == species]
    protein_names = sorted(list(set([i[4] for i in db_proteins])))
    db_fields = [i for i in db_fields if i[3] in protein_names]
    protein_number_before = (len(list(set([i[3] for i in db_fields]))))
    print("{0} were performed in {1}.\n".format(len(db_fields), species))
    db_fields = [i for i in db_fields if i[2] != ""]
    protein_number_after = (len(list(set([i[3] for i in db_fields]))))
    db_fields = [[i[3], "RBPDB", i[0], i[1], i[2]] for i in db_fields]
    print("After removing experiments with no reported motif, {0} proteins remain of the initial {1}.\n".format(protein_number_after, protein_number_before))

    bases = np.array(["A", "C", "G", "U"])
    db_pwm_list = rw.read_many_fields(RBPDB_PWMs, "\t")

    for i in db_pwm_list:
        if i[1] in protein_names:
            current_file_name = "{0}/{1}.pwm".format(pwm_dir, i[0])
            current_PWM = rw.read_many_fields(current_file_name, delimiter = " ")
            for j in range(len(current_PWM)):
                current_PWM[j] = [float(k) for k in current_PWM[j] if k != ""]
            consensus = nc.consensus_from_PWM(current_PWM, bases, 0)
            PMID = i[0].split("_")
            PMID = PMID[1]
            new_record = [i[1], "RBPDB_PWM", PMID, "SELEX", consensus]
            db_fields.append(new_record)

    protein_number_after = (len(list(set([i[0] for i in db_fields]))))
    print("After adding additional sequences from SELEX PWMs (RBPDB), there are {0} proteins.\n".format(protein_number_after))

    if species == "Mus musculus":
        RBPmap_proteins = rw.read_many_fields("RBP/RBPmap_proteins.csv", ",")
        RBPmap_proteins = list_to_dict(RBPmap_proteins, 0, 1)
        RNAc_source = [i for i in RBPmap_proteins if "23846655" in RBPmap_proteins[i]]
    else:
        RNAc_source = []

    for file_name in os.listdir(RBPmap_PSSMs):
        #RBPmap and SFmap don't distinguish between human and mouse motifs
        if "human" in file_name:
            file_name_split = file_name.split("_")
            protein_name = file_name_split[0]
            if protein_name not in RNAc_source:
                initial_pssm = rw.read_many_fields(os.path.join(RBPmap_PSSMs, file_name), delimiter = "\t")
                current_pssm = initial_pssm[1:]
                current_pssm = [i[1:] for i in current_pssm]
                for i in range(len(current_pssm)):
                    current_pssm[i] = [float(j) for j in current_pssm[i]]
                consensus = nc.consensus_from_PWM(current_pssm, bases, 0.25, transform = True)
                protein_name = list(protein_name)
                if protein_name[:4] == ["S", "R", "S", "F"]:
                    protein_name[:4] = ["S", "F", "R", "S"]
                protein_name = "".join(protein_name)
                new_record = [protein_name, "RBPmap_PWM", "NULL", "various", consensus]
                db_fields.append(new_record)

    protein_number_after = (len(list(set([i[0] for i in db_fields]))))
    print("After adding additional sequences from RBPmap PSSMs, there are {0} proteins.\n".format(protein_number_after))

    SFmap_data = rw.read_many_fields(SFmap_proteins, delimiter = ",")

    for i in SFmap_data:
        if "," in i[1]:
            temp_split = i[1].split(", ")
            temp_split = [j.upper() for j in temp_split]
            i[1] = ";".join(temp_split)
        else:
            i[1] = i[1].upper()
        new_record = [i[0], "SFmap", "NULL", "various", i[1]]
        db_fields.append(new_record)

    protein_number_after = (len(list(set([i[0] for i in db_fields]))))
    print("After adding motifs from SFmap, there are {0} proteins.\n".format(protein_number_after))

    RNAc = rw.read_many_fields(RNAcompete_information, delimiter = "\t")
    RNAc = [i for i in RNAc[1:] if i]
    if species == "Homo sapiens":
        RNAc = [i for i in RNAc if i[3] != "." and i[8] == "D"]
    if species == "Mus musculus":
        RNAc = [i for i in RNAc if i[3] != "."]

    PSSM_folder = RNAcompete_PWMs
    for record in RNAc:
        motif_name = record[3]
        initial_pssm = rw.read_many_fields(os.path.join(PSSM_folder, "{0}.txt".format(motif_name)), delimiter = "\t")
        if initial_pssm == []:
            if record[19] == "21036867":#RBPDB paper
                pass
            else:
                print(record)
        else:    
            current_pssm = initial_pssm[1:]
            current_pssm = [i[1:] for i in current_pssm]
            for i in range(len(current_pssm)):
                current_pssm[i] = [float(j) for j in current_pssm[i]]
            consensus = nc.consensus_from_PWM(current_pssm, bases, 0.25, transform = True)
            protein_name = record[6]
            new_record = [protein_name, "CIS-BP_RNA_PWM", record[19], record[14], consensus] 
            db_fields.append(new_record)

    protein_number_after = (len(list(set([i[0] for i in db_fields]))))
    print("After adding motifs from CIS-BP RNA, there are {0} proteins.\n".format(protein_number_after))

    to_delete = []
    for pos, i in enumerate(db_fields):
        if ";" in i[4]:
            if "; " in i[4]:
                temp_split = i[4].split("; ")
            else:
                temp_split = i[4].split(";")
            temp_split = [((j.upper()).lstrip("N")).rstrip("N") for j in temp_split]
            temp_split = [j for j in temp_split if len(j) <= upper_threshold and len(j) >= lower_threshold and "(" not in j]
            if temp_split:
                db_fields[pos][4] = temp_split[0]
                for j in temp_split[1:]:
                    db_fields.append([i[0], i[1], i[2], i[3], j])
            else:
                to_delete.append(pos)
        else:
            i[4] = (((i[4]).upper()).rstrip("N")).lstrip("N")
            if len(i[4]) > upper_threshold or len(i[4]) < lower_threshold or "(" in i[4]:
                to_delete.append(pos)
            else:
                db_fields[pos][4] = i[4]

    db_fields = [i for pos, i in enumerate(db_fields) if pos not in to_delete]

    protein_number_after = (len(list(set([i[0] for i in db_fields]))))
    print("After only keeping motifs of length {0}-{1} bp, {2} proteins remain.\n".format(lower_threshold, upper_threshold, protein_number_after))

    protein_names = list(set([i[0] for i in db_fields]))

    if species == "Mus musculus":
        protein_names_file = "RBP/RBP_names_for_checking.txt"
        with open(protein_names_file, "w") as file:
            for name in protein_names:
                file.write("{0}\n".format(name))
        MGI_file = "RBP/MGI_correspondances.txt"
        MGI = rw.read_many_fields(MGI_file, "\t")
        MGI_names_all = [i[0] for i in MGI[1:]]
        found = [i[0] for i in MGI if i[0] == i[3]]
        MGI = {i[0]: i[3] for i in MGI[1:] if i[0] not in found}

    to_delete = []
    for pos, i in enumerate(db_fields):
        if species == "Mus musculus":
            db_fields[pos][0] = "".join([db_fields[pos][0][0].upper(), db_fields[pos][0][1:].lower()])
            #will get rid of Hnrnpcl1, which didn't return anything in the MGI search.
            if db_fields[pos][0] not in MGI_names_all:
                to_delete.append(pos)
            else:
                if db_fields[pos][0] not in found:
                    db_fields[pos][0] = MGI[db_fields[pos][0]]
        elif species == "Homo sapiens":
            if i[0] == "A2BP1" or i[0] == "FOX1":
                db_fields[pos][0] = "RBFOX1"
            elif i[0] == "SFRS13A":
                db_fields[pos][0] = "SRSF10"
            elif i[0][:6] == "BRUNOL":
                db_fields[pos][0] = "CELF{0}".format(i[0][-1])
            elif i[0] == "CUGBP":
                db_fields[pos][0] = "CELF1"
            elif i[0] == "Fusip1":
                db_fields[pos][0] = "SRSF10"
            elif i[0][:4] == "SFRS":
                db_fields[pos][0] = "SRSF{0}".format(i[0][4:])
            elif i[0] == "HuR":
                db_fields[pos][0] = "ELAVL1"
            elif i[0] == "MBNL":
                db_fields[pos][0] = "MBNL1"
            elif i[0] == "PTB":
                db_fields[pos][0] = "PTBP1"
            elif i[0] == "QK1":
                db_fields[pos][0] = "QKI"
            elif i[0] == "RBM9":
                db_fields[pos][0] = "RBFOX2"
            elif i[0] == "STAR-PAP":
                db_fields[pos][0] = "TUT1"
            elif i[0] == "YB-1":
                db_fields[pos][0] = "YBX1"
            elif i[0] == "hnRNPK":
                db_fields[pos][0] = "HNRNPK"
            elif i[0] == "hnRNPLL" or i[0] == "HNRPLL":
                db_fields[pos][0] = "HNRNPLL"

    db_fields = [i for pos, i in enumerate(db_fields) if pos not in to_delete]

    protein_names = list(set([i[0] for i in db_fields]))

    protein_number_after = (len(list(set([i[0] for i in db_fields]))))
    print("After cleaning up protein IDs, {0} proteins remain.\n".format(protein_number_after))
            
    protein_dict = {}
    for i in db_fields:
        if i[0] not in protein_dict.keys():
            protein_dict[i[0]] = [i]
        else:
            protein_dict[i[0]].append(i)

    if species == "Homo sapeins":
        del protein_dict["PPIE"]
        del protein_dict["MIR1236"]
        del protein_dict["PABPC4"]
        print("After removing PPIE, PABPC4 and MIR1236, {0} proteins remain.\n".format(len(protein_dict)))
    elif species == "Mus musculus":
        del protein_dict["Pabpc4"]
        print("After removing Pabpc4, {0} proteins remain.\n".format(len(protein_dict)))

    for i in protein_dict:
        if i == "ELAVL1":
            protein_dict[i].append(['ELAVL1', 'synthetic', 'synthetic', 'synthetic', 'UUWGDUU'])
        elif i == "ELAVL2":
            protein_dict[i].append(['ELAVL2', 'synthetic', 'synthetic', 'synthetic', 'RWUUYAUUUWR'])
        protein_dict[i] = sorted(protein_dict[i], key = lambda x:x[4])
        current_motifs = [j[4] for j in protein_dict[i]]
        to_delete = []
        for j in range(1, len(current_motifs)):
            if current_motifs[j] == current_motifs[j-1]:
                for k in range(1, 4):
                    protein_dict[i][j][k] = ",".join([protein_dict[i][j][k], protein_dict[i][j - 1][k]])
                to_delete.append(j - 1)
        protein_dict[i] = [protein_dict[i][j] for j in range(len(protein_dict[i])) if j not in to_delete]

    for i in protein_dict:
        protein_dict[i] = [[j[0], j[4], j[1], j[2], j[3]] for j in protein_dict[i]] 

    print("\n")
    print("Writing motifs to {0}.\n".format(final_motifs_file_name))

    motif_numbers = []
    with open(final_motifs_file_name, "w") as final_motifs_file:
        for i in sorted(list(protein_dict.keys())):
            final_motifs_file.write(">{0}\n".format(i))
            current_motifs = [j[1] for j in protein_dict[i]]
            DNA_motifs = [nc.DNA_RNA_conversion(j) for j in current_motifs]
            unravelled_motifs = [nc.unravel_consensus(j) for j in DNA_motifs]
            unravelled_motifs = flatten(unravelled_motifs)
            unravelled_motifs = list(set(unravelled_motifs))
            print("Writing {0} motifs for {1}.".format(len(unravelled_motifs), i))
            motif_numbers.append(len(unravelled_motifs))
            unravelled_motifs = "|".join(unravelled_motifs)
            final_motifs_file.write("{0}\n".format(unravelled_motifs))

    plt.figure(1)
    plotting.histogram(motif_numbers, 50, x_lab = "Motif number", y_lab = "Frequency", title = None)
    plotting.save_and_show([10, 10], 100, plot_name)

if __name__ == "__main__":
    main()





