from django.shortcuts import render
import requests

import primer3
from Bio import Phylo
import pylab
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import AlignIO
from Bio import Entrez
import re
from Bio import SeqIO
def home(request):
    return render(request, 'index.html')
def snp(request):
    rs = request.GET.get('rs' ,'rs141033578')
    clinical_ann_metadata = open('staticfiles/clinical_ann_metadata.tsv', 'r').readlines()
    clinical_ann = open('staticfiles/clinical_ann.tsv', 'r').readlines()
    var_drug_ann = open('staticfiles/var_drug_ann.tsv','r').readlines()
    study_parameters = open('staticfiles/study_parameters.tsv','r').readlines()


    # Clinical annotation metadata file...!
    clinical_annotation_ID = []
    Location = []
    Gene = []
    level_of_evidence = []
    Clinical_Annotation_Types = []
    Genotype_Phenotype_IDs = []
    annotation_text = []
    Variant_Annotations_IDs = []
    Variant_Annotations = []
    Evidence_Count = []
    Related_Chemicals = []
    Related_Diseases = []
    Race = []
    Chromosome = []
    pMID = []
    other_diseases=[]
    for linecam in clinical_ann_metadata:
        linecam = re.split('\t', linecam)
        # print(line[1])
        if rs in linecam[1]:
            try:
                clinical_annotation_ID.append(linecam[0])
                Location.append(linecam[1])
                linecam[2] = re.split('\s', linecam[2])
                genename=linecam[2][0]
                Gene.append(linecam[2][0])
                level_of_evidence.append(linecam[3])
                Clinical_Annotation_Types.append(linecam[4])
                Genotype_Phenotype_IDs.append(linecam[5])
                annotation_text.append("%s\n"%linecam[6])
                Variant_Annotations_IDs.append(linecam[7])
                Variant_Annotations.append(linecam[8])
                pMID.append(linecam[9])
                Evidence_Count.append(linecam[10])
                linecam[11] = re.split('\s', linecam[11])
                Related_Chemicals.append(linecam[11][0])
                Related_Diseases.append(linecam[12])
                disaseHTML=linecam[12]
                disaseHTML=re.split('\s',disaseHTML)
                print(disaseHTML[0])
                print(disaseHTML[1])
                disease_name="%s %s"%(disaseHTML[0].replace('"',''),disaseHTML[1].replace('"',''))
                print(disease_name)
                other_diseases.append(disease_name)
                print('*****************************')
                Race.append(linecam[13])
                Chromosome.append(linecam[14])
            except:
                pass



    linecam5 = re.split(',', linecam[5])
    linenum = 0
    # clinical_ann
    genophenoID = []
    clinicalpheno = []
    genotype = []
    for lineca in clinical_ann:
        lineca = re.split('\t', lineca)
        try:
            if linecam5[int(linenum)].replace('"', '') in lineca:
                genophenoID.append(linecam5[int(linenum)])
                clinicalpheno.append(lineca[2])
                genotype.append(lineca[1])
                linenum += 1
        except:
            pass
    # var_drug_ann file....
    Annotation_ID = []
    Variant = []
    Genevda = []
    Chemical = []
    PMIDva = []
    Phenotype_Category = []
    Significance = []
    Notes = []
    Sentence = []
    StudyParameters = []
    Alleles = []
    Chromosomevda = []
    for linevda in var_drug_ann:
        linevda = re.split('\t', linevda)
        if linecam[2] in linevda[2]:
            Annotation_ID.append(linevda[0])
            Variant.append(linevda[1])
            Genevda.append(linevda[2])
            Chemical.append(linevda[3])
            PMIDva.append(linevda[4])
            Phenotype_Category.append(linevda[5])
            Significance.append(linevda[6])
            Notes.append(linevda[7])
            Sentence.append(linevda[8])
            StudyParameters.append(linevda[9])
            Alleles.append(linevda[10])
            Chromosomevda.append(linevda[11])

    # Study_parameters
    Study_Parameters_ID = []
    Study_Type = []
    Study_Cases = []
    Study_Controls = []
    Characteristics = []
    Characteristics_Type = []
    Frequency_In_Cases = []
    Allele_Of_Frequency_In_Cases = []
    Frequency_In_Controls = []
    Allele_Of_Frequency_In_Controls = []
    P_Value_Operator = []
    P_Value = []
    Ratio_Stat_Type = []
    Ratio_Stat = []
    Confidence_Interval_Start = []
    Confidence_Interval_Stop = []
    Race_s = []
    num=[]

    for i in range(100):
        num.append(i)

    for linesp in study_parameters:
        linesp = re.split('\t', linesp)
        if linevda[9] in linesp[0]:
            Study_Parameters_ID.append(linesp[0])
            Study_Type.append(linesp[1])
            Study_Cases.append(linesp[2])
            Study_Controls.append(linesp[3])
            Characteristics.append(linesp[4])
            Characteristics_Type.append(linesp[5])
            Frequency_In_Cases.append(linesp[6])
            Allele_Of_Frequency_In_Cases.append(linesp[7])
            Frequency_In_Controls.append(linesp[8])
            Allele_Of_Frequency_In_Controls.append(linesp[9])
            P_Value_Operator.append(linesp[10])
            P_Value.append(linesp[11])
            Ratio_Stat_Type.append(linesp[12])
            Ratio_Stat.append(linesp[13])
            Confidence_Interval_Start.append(linesp[14])
            Confidence_Interval_Stop.append(linesp[15])
            Race_s.append(linesp[16])
    annotation_text = [i.replace('"', '') for i in annotation_text]
    #
    # alllists = zip(clinical_annotation_ID, Location, Gene, level_of_evidence, Clinical_Annotation_Types,
    #                Genotype_Phenotype_IDs, annotation_text, Variant_Annotations_IDs[0:2], Variant_Annotations[:1],
    #                Evidence_Count, Related_Chemicals, Related_Diseases, Race, Chromosome, pMID)
    #
    # alllists = zip(clinical_annotation_ID, Location, Gene, level_of_evidence, Clinical_Annotation_Types, Genotype_Phenotype_IDs, annotation_text, Variant_Annotations_IDs, Variant_Annotations, Evidence_Count, Related_Chemicals, Related_Diseases, Race, Chromosome, pMID, Study_Type, Study_Cases, Study_Controls, Characteristics, Characteristics_Type, Frequency_In_Cases, Allele_Of_Frequency_In_Cases, Frequency_In_Controls, Allele_Of_Frequency_In_Controls, P_Value_Operator, P_Value, Ratio_Stat_Type, Ratio_Stat, Confidence_Interval_Start, Confidence_Interval_Stop, Race_s)



    # print(clinical_annotation_ID)

    reactom_id_option=[]
    reactom_id_info=[]

    geneid = open('staticfiles/clinical_ann_metadata.tsv','r').readlines()
    genencbi = open('staticfiles/genes.tsv','r').readlines()
    file = open('staticfiles/Ensembl_Pathway.tsv','r').readlines()


    for line in geneid:
        line = re.split('\s', line)
        if rs == line[1]:
            gene_name = line[2]
            print(gene_name)
            for id in genencbi:
                id = re.split('\t', id)
                # print(id[5])
                if gene_name == id[5]:
                    print(id[1] + "**********")
                    print(id[0] + '----------')
                    ncbiid = id[1]
                    print(ncbiid)
                    try:
                        ensmbleid = id[3]
                    except UnboundLocalError:
                        ensmbleid = ''
                    print(ensmbleid)
                    handle = Entrez.efetch(db="nuccore", id='%s' % ncbiid, rettype="gb", retmode="text")
                    record = SeqIO.parse(handle, "genbank")
                    for rec in record:
                        print('')
                        print('')
                        keyw = rec.annotations["keywords"]
                        tax = rec.annotations["taxonomy"]

            for line in file:
                line = re.split('\t', line)
                if ensmbleid in line:
                    # print(line)
                    if line[3].startswith(" Defective"):
                        reactomid = line[1]
                        print(reactomid)
                        # print(line)
                        break
                    elif gene_name == line[3]:
                        reactomid = line[1]
                        print(reactomid)
                        # print(line)
                        break
                    elif gene_name in line[5]:
                        reactomid = line[1]
                        print(reactomid)
                        # print(line)
                        break

                    else:
                        reactomid = ''

            for line in file:
                line = re.split('\t', line)
                if ensmbleid in line:
                    # print(line)
                    if gene_name in line[3]:
                        reactomid_to_list = line[1]
                        reactomid_info = line[3]
                        reactom_id_option.append(reactomid_to_list)
                        reactom_id_info.append(reactomid_info)
                        print(line)
                    else:
                        reactomid_to_list = ''


            url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id='
            out = open('out.txt', 'w')
            response = requests.get(url+ncbiid)
            out.write(response.text)
            out.close()
            out = open('out.txt', 'r').readlines()
            for line in out:
                print(line)
            for line in out:
                if '<Name>' in line:
                    name = line.strip().replace('<Name>', '').replace('</Name>', '')
                if '<NomenclatureName>' in line:
                    nomenclature = line.strip().replace('<NomenclatureName>', '').replace('</NomenclatureName>', '')
                if '<Description>' in line:
                    description = line.strip().replace('<Description>', '').replace('</Description>', '')
                if '<MapLocation>' in line:
                    maplocation = line.strip().replace('<MapLocation>', '').replace('</MapLocation>', '')
                if '<OtherAliases>' in line:
                    otherAliases = line.strip().replace('<OtherAliases>', '').replace('</OtherAliases>', '')
                if '<ExonCount>' in line:
                    exoncount = line.strip().replace('<ExonCount>', '').replace('</ExonCount>', '')
                if '<GeneWeight>' in line:
                    geneweight = line.strip().replace('<GeneWeight>', '').replace('</GeneWeight>', '')
                if '<Summary>' in line:
                    summary = line.strip().replace('<Summary>', '').replace('</Summary>', '')
                if '<ScientificName>' in line:
                    sciname = line.strip().replace('<ScientificName>', '').replace('</ScientificName>', '')
                if '<OtherDesignations>' in line:
                    designation = line.strip().replace('<OtherDesignations>', '').replace('</OtherDesignations>', '')
                if '<ChrLoc>' in line:
                    chrloc = line.strip().replace('<ChrLoc>', '').replace('</ChrLoc>', '')
                if '<ChrAccVer>' in line:
                    chracc = line.strip().replace('<ChrAccVer>', '').replace('</ChrAccVer>', '')
                if '<ChrStart>' in line:
                    chrstart = line.strip().replace('<ChrStart>', '').replace('</ChrStart>', '')
                if '<ChrStop>' in line:
                    chrend = line.strip().replace('<ChrStop>', '').replace('</ChrStop>', '')

            # 'apiname': name, 'description': description, 'maplocation': maplocation,
            # 'otherAliases': otherAliases, 'exoncount': exoncount, 'geneweight': geneweight,
            # 'summary': summary, 'sciname': sciname, 'designation': designation, 'chrloc': chrloc, 'chracc': chracc, 'chrstart': chrstart, 'chrend': chrend



    reactome_list=zip(reactom_id_option,reactom_id_info)
    clinical= zip(clinical_annotation_ID, Location, Gene, level_of_evidence, Clinical_Annotation_Types, Genotype_Phenotype_IDs, annotation_text, Variant_Annotations_IDs, Variant_Annotations, Evidence_Count, Related_Chemicals, Related_Diseases, Race, Chromosome, pMID,num)
    study= zip(Study_Type, Study_Cases, Study_Controls, Characteristics, Characteristics_Type, Frequency_In_Cases, Allele_Of_Frequency_In_Cases, Frequency_In_Controls, Allele_Of_Frequency_In_Controls, P_Value_Operator, P_Value, Ratio_Stat_Type, Ratio_Stat, Confidence_Interval_Start, Confidence_Interval_Stop, Race_s)
    vardrugann= zip(Annotation_ID,Variant,Genevda,Chemical,PMIDva,Phenotype_Category,Significance,Notes,Sentence,StudyParameters,Alleles,Chromosomevda)





    all_reactom=zip(reactom_id_option,reactom_id_info)
    return render(request, 'snp(gencardform).html', context={'clinical':clinical,
                                                             'study':study,'snpname':rs,'desc':rec.description,'org':rec.annotations["organism"],
                                                             'date':rec.annotations["date"],'keyw':keyw,'src':rec.annotations["source"],'tax':tax,'id':rec.id,
                                                             'genename':genename,'ncbiid':ncbiid,'seq':rec.seq[:50]+'....','disease':disaseHTML,'ensid':ensmbleid,
                                                             'reactomid':reactomid,'all_reactom':all_reactom,'reactome_list':reactome_list,'rs':rs,
                                                             'gene_name':gene_name,'other_diseases':other_diseases,'apiname': name, 'description': description, 'maplocation': maplocation,
                                                             'otherAliases': otherAliases, 'exoncount': exoncount,
                                                             'geneweight': geneweight,'summary': summary, 'sciname': sciname, 'designation': designation,
                                                             'chrloc': chrloc, 'chracc': chracc, 'chrstart': chrstart, 'chrend': chrend,'nomenclature':nomenclature})
def primerdesighn(request, primer3=None):
    snpid = request.GET.get('rs' ,'rs141033578')
    SEQUENCE_INCLUDED_REGION_START = request.GET.get('SEQUENCE_INCLUDED_REGION_START', 30)
    PRIMER_OPT_SIZE = request.GET.get('PRIMER_OPT_SIZE', 20)
    PRIMER_PICK_INTERNAL_OLIGO = request.GET.get('PRIMER_PICK_INTERNAL_OLIGO', 1)
    PRIMER_INTERNAL_MAX_SELF_END = request.GET.get('PRIMER_INTERNAL_MAX_SELF_END', 8)
    PRIMER_MIN_SIZE = request.GET.get('PRIMER_MIN_SIZE', '18')
    PRIMER_MAX_SIZE = request.GET.get('PRIMER_MAX_SIZE', 25)
    PRIMER_OPT_TM = request.GET.get('PRIMER_OPT_TM', 60)
    PRIMER_MIN_TM = request.GET.get('PRIMER_MIN_TM', 57)
    PRIMER_MAX_TM = request.GET.get('PRIMER_MAX_TM', 63)
    PRIMER_MIN_GC = request.GET.get('PRIMER_MIN_GC', 20)
    PRIMER_MAX_GC = request.GET.get('PRIMER_MAX_GC', 80)
    PRIMER_MAX_POLY_X = request.GET.get('PRIMER_MAX_POLY_X', 100)
    PRIMER_INTERNAL_MAX_POLY_X = request.GET.get('PRIMER_INTERNAL_MAX_POLY_X', 100)
    PRIMER_SALT_MONOVALENT = request.GET.get('PRIMER_SALT_MONOVALENT', 50)
    PRIMER_DNA_CONC = request.GET.get('PRIMER_DNA_CONC', 50)
    PRIMER_MAX_NS_ACCEPTED = request.GET.get('PRIMER_MAX_NS_ACCEPTED', 0)
    PRIMER_MAX_SELF_ANY = request.GET.get('PRIMER_MAX_SELF_ANY', 12)
    PRIMER_MAX_SELF_END = request.GET.get('PRIMER_MAX_SELF_END', 8)
    PRIMER_PAIR_MAX_COMPL_ANY = request.GET.get('PRIMER_PAIR_MAX_COMPL_ANY', 12)
    PRIMER_PAIR_MAX_COMPL_END = request.GET.get('PRIMER_PAIR_MAX_COMPL_END', 8)
    PRIMER_PRODUCT_SIZE_RANGE_START = request.GET.get('PRIMER_PRODUCT_SIZE_RANGE_START', 27)
    PRIMER_PRODUCT_SIZE_RANGE_END = request.GET.get('PRIMER_PRODUCT_SIZE_RANGE_END', 100000000000)

    info = []
    results = []
    seqs = []
    Entrez.email = "A.N.Other@example.com"  # Always tell NCBI who you are
    entrez_id = "%s" % snpid
    handle = Entrez.efetch(db="snp", id=snpid, rettype="fasta", retmode="fasta")
    record = SeqIO.read(handle, "fasta")
    SeqIO.write(record, "" + entrez_id + ".fasta", "fasta")
    record = SeqIO.read("" + entrez_id + ".fasta", "fasta")
    # print("%s %i" % (record.seq, int(len(record))))
    sequence = record.seq
    file = open('seq.txt', 'w')
    file.write(str(sequence))
    file.close()
    file = open('seq.txt', 'r').readlines()
    sequence = ''
    for line in file:
        line = line.replace('R', '').replace('B','').replace('H','').replace('K','')
        # print(line)
        sequence = line
    print(sequence)
    seqlen = len(record)-10
    SEQUENCE_INCLUDED_REGION_END = request.GET.get('SEQUENCE_INCLUDED_REGION_END', seqlen-100)
    seqid = record.id
    primer = primer3.bindings.designPrimers(
        {
            'SEQUENCE_ID': seqid,
            'SEQUENCE_TEMPLATE': '%s' % sequence,
            'SEQUENCE_INCLUDED_REGION': [0,seqlen-10]
        },
        {
            'PRIMER_OPT_SIZE': int(int(PRIMER_OPT_SIZE)),
            'PRIMER_PICK_INTERNAL_OLIGO': int(int(PRIMER_PICK_INTERNAL_OLIGO)),
            'PRIMER_INTERNAL_MAX_SELF_END': int(PRIMER_INTERNAL_MAX_SELF_END),
            'PRIMER_MIN_SIZE': int(PRIMER_MIN_SIZE),
            'PRIMER_MAX_SIZE': int(PRIMER_MAX_SIZE),
            'PRIMER_OPT_TM': int(PRIMER_OPT_TM),
            'PRIMER_MIN_TM': int(PRIMER_MIN_TM),
            'PRIMER_MAX_TM': int(PRIMER_MAX_TM),
            'PRIMER_MIN_GC': int(PRIMER_MIN_GC),
            'PRIMER_MAX_GC': int(PRIMER_MAX_GC),
            'PRIMER_MAX_POLY_X': int(PRIMER_MAX_POLY_X),
            'PRIMER_INTERNAL_MAX_POLY_X': int(PRIMER_INTERNAL_MAX_POLY_X),
            'PRIMER_SALT_MONOVALENT': int(PRIMER_SALT_MONOVALENT),
            'PRIMER_DNA_CONC': int(PRIMER_DNA_CONC),
            'PRIMER_MAX_NS_ACCEPTED': int(PRIMER_MAX_NS_ACCEPTED),
            'PRIMER_MAX_SELF_ANY': int(PRIMER_MAX_SELF_ANY),
            'PRIMER_MAX_SELF_END': int(PRIMER_MAX_SELF_END),
            'PRIMER_PAIR_MAX_COMPL_ANY': int(PRIMER_PAIR_MAX_COMPL_ANY),
            'PRIMER_PAIR_MAX_COMPL_END': int(PRIMER_PAIR_MAX_COMPL_END),
            'PRIMER_PRODUCT_SIZE_RANGE': [[int(PRIMER_PRODUCT_SIZE_RANGE_START), int(PRIMER_PRODUCT_SIZE_RANGE_END)]],
        })

    # for i in primer:
    #     print(i, primer[i])
    primerinfo0 = []
    primerresult0 = []
    primerinfo1 = []
    primerresult1 = []
    primerinfo2 = []
    primerresult2 = []
    primerinfo3 = []
    primerresult3 = []
    primerinfo4 = []
    primerresult4 = []

    for key in sorted(primer):
        # key = "%s: %s" % (key, primer[key])
        if '_0_' in key:
            primerresult0.append(primer[key])
            if key in key:
                key=key.replace('0','').replace('_',' ').replace('th','TH').replace('tm','TM')
                key=key.lower().capitalize()
                primerinfo0.append(key)


        elif '_1_' in key:
            primerinfo1.append(key)
            primerresult1.append(primer[key])

        elif '_2_' in key:
            primerinfo2.append(key)
            primerresult2.append(primer[key])

        elif '_3_' in key:
            primerinfo3.append(key)
            primerresult3.append(primer[key])


        elif '_4_' in key:
            primerinfo4.append(key)
            primerresult4.append(primer[key])
    # alldata = zip(info, results)
    # print(alldata)

    # print("%s\n%s\n%s\n%s\n%s\n"%(primerinfo0,primerinfo1,primerinfo2,primerinfo3,primerinfo4))


    alldata = zip(primerinfo0, primerresult0, primerinfo1, primerresult1, primerinfo2, primerresult2, primerinfo3, primerresult3,primerinfo4, primerresult4)
    print(alldata)
    primer1= zip(primerinfo0,primerresult0)
    primer2= zip(primerinfo1,primerresult1)
    primer3= zip(primerinfo2,primerresult2)
    primer4= zip(primerinfo3,primerresult3)
    primer5= zip(primerinfo4,primerresult4)
    return render(request, 'pick_primer.html',context={'alldata':alldata,'length':seqlen})




def drug(request):
    import re
    PharmGKB_Accession_Id = []
    Name = []
    Generic_Names = []
    Trade_Names = []
    Brand_Mixtures = []
    Type = []
    Cross_references = []
    SMILES = []
    InChI = []
    Dosing_Guideline = []
    External_Vocabulary = []
    Clinical_Annotation_Count = []
    Variant_Annotation_Count = []
    Pathway_Count = []
    VIP_Count = []
    Dosing_Guideline_Sources = []
    Top_Clinical_Annotation_Level = []
    Top_FDA_Label_Testing_Level = []
    Top_Any_Drug_Label_Testing_Level = []
    Label_Has_Dosing_Info = []
    Has_Rx_Annotation = []
    num=[]
    inpo =  request.GET.get('drug' ,'Lipid')

    file = open('staticfiles/drugs.tsv',
                'r').readlines()
    for line in file:
        line = re.split('\t', line.strip())
        if inpo in line[1]:
            try:
                PharmGKB_Accession_Id.append(line[0])
                Name.append(line[1])
                Generic_Names.append(line[2])
                Trade_Names.append(line[3])
                Brand_Mixtures.append(line[4])
                Type.append(line[5])
                Cross_references.append(line[6])
                SMILES.append(line[7])
                InChI.append(line[8])
                Dosing_Guideline.append(line[9])
                External_Vocabulary.append(line[10])
                Clinical_Annotation_Count.append(line[11])
                Variant_Annotation_Count.append(line[12])
                Pathway_Count.append(line[13])
                VIP_Count.append(line[14])
                Dosing_Guideline_Sources.append(line[15])
                Top_Clinical_Annotation_Level.append(line[16])
                Top_FDA_Label_Testing_Level.append(line[17])
                Top_Any_Drug_Label_Testing_Level.append(line[18])
                Label_Has_Dosing_Info.append(line[19])
                Has_Rx_Annotation.append(line[20])
            except:
                continue

    for val in range(100):
        num.append(str(val))

    drug = zip(PharmGKB_Accession_Id, Name, Generic_Names, Trade_Names, Brand_Mixtures, Type, Cross_references, SMILES,
               InChI, Dosing_Guideline, External_Vocabulary, Clinical_Annotation_Count, Variant_Annotation_Count,
               Pathway_Count, VIP_Count, Dosing_Guideline_Sources, Top_Clinical_Annotation_Level,
               Top_FDA_Label_Testing_Level, Top_Any_Drug_Label_Testing_Level, Label_Has_Dosing_Info, Has_Rx_Annotation , num)
    print(Name)
    return render(request, 'drug.html', context={'drug':drug})
def chemical(request):
    import re
    PharmGKB_Accession_Id = []
    Name = []
    Generic_Names = []
    Trade_Names = []
    Brand_Mixtures = []
    Type = []
    Cross_references = []
    SMILES = []
    InChI = []
    Dosing_Guideline = []
    External_Vocabulary = []
    Clinical_Annotation_Count = []
    Variant_Annotation_Count = []
    Pathway_Count = []
    VIP_Count = []
    Dosing_Guideline_Sources = []
    Top_Clinical_Annotation_Level = []
    Top_FDA_Label_Testing_Level = []
    Top_Any_Drug_Label_Testing_Level = []
    Label_Has_Dosing_Info = []
    Has_Rx_Annotation = []
    num=[]
    inpo =  request.GET.get('chemical' ,'ivacaftor')

    file = open('staticfiles/chemicals.tsv',
                'r').readlines()
    for line in file:
        line = re.split('\t', line.strip())
        if inpo in line[1]:
            try:
                PharmGKB_Accession_Id.append(line[0])
                Name.append(line[1])
                Generic_Names.append(line[2])
                Trade_Names.append(line[3])
                Brand_Mixtures.append(line[4])
                Type.append(line[5])
                Cross_references.append(line[6])
                SMILES.append(line[7])
                InChI.append(line[8])
                Dosing_Guideline.append(line[9])
                External_Vocabulary.append(line[10])
                Clinical_Annotation_Count.append(line[11])
                Variant_Annotation_Count.append(line[12])
                Pathway_Count.append(line[13])
                VIP_Count.append(line[14])
                Dosing_Guideline_Sources.append(line[15])
                Top_Clinical_Annotation_Level.append(line[16])
                Top_FDA_Label_Testing_Level.append(line[17])
                Top_Any_Drug_Label_Testing_Level.append(line[18])
                Label_Has_Dosing_Info.append(line[19])
                Has_Rx_Annotation.append(line[20])
            except:
                pass

    for val in range(100):
        num.append(str(val))

    drug = zip(PharmGKB_Accession_Id, Name, Generic_Names, Trade_Names, Brand_Mixtures, Type, Cross_references, SMILES,
               InChI, Dosing_Guideline, External_Vocabulary, Clinical_Annotation_Count, Variant_Annotation_Count,
               Pathway_Count, VIP_Count, Dosing_Guideline_Sources, Top_Clinical_Annotation_Level,
               Top_FDA_Label_Testing_Level, Top_Any_Drug_Label_Testing_Level, Label_Has_Dosing_Info, Has_Rx_Annotation , num)
    print(Name)
    # print(num)
    return render(request, 'chemical.html', context={'drug':drug})


def genes(request):
    import re
    import requests
    import re

    inpo = request.GET.get('genes', 'CFTR')
    PharmGKB_Accession_Id = []
    NCBI_Gene_ID = []
    HGNC_ID = []
    Ensembl_Id = []
    Name = []
    Symbol = []
    Alternate_Names = []
    Alternate_Symbols = []
    Is_VIP = []
    Has_Variant_Annotation = []
    Cross_references = []
    Has_CPIC_Dosing_Guideline = []
    Chromosome = []
    Chromosomal_Start_GRCh37_p13 = []
    Chromosomal_Stop_GRCh37_p13 = []
    Chromosomal_Start_GRCh38_p7 = []
    Chromosomal_Stop_GRCh38_p7 = []

    file = open('staticfiles/genes.tsv',
                'r').readlines()
    for line in file:
        line = re.split('\t', line)
        if inpo == line[5]:
            PharmGKB_Accession_Id.append(line[0])
            NCBI_Gene_ID.append(line[1])
            ncbiid = line[1]
            HGNC_ID.append(line[2])
            Ensembl_Id.append(line[3])
            ensmbleid = line[3]
            Name.append(line[4])
            Symbol.append(line[5])
            namegene=line[5]
            Alternate_Names.append(line[6])
            Alternate_Symbols.append(line[7])
            Is_VIP.append(line[8])
            Has_Variant_Annotation.append(line[9])
            Cross_references.append(line[10])
            Has_CPIC_Dosing_Guideline.append(line[11])
            Chromosome.append(line[12])
            Chromosomal_Start_GRCh37_p13.append(line[13])
            Chromosomal_Stop_GRCh37_p13.append(line[14])
            Chromosomal_Start_GRCh38_p7.append(line[15])
            Chromosomal_Stop_GRCh38_p7.append(line[16])
            print(ncbiid)

            file = open('staticfiles/Ensembl_Pathway.tsv',
                        'r').readlines()

            for line in file:
                line = re.split('\t', line)
                if ensmbleid in line:
                    # print(line)
                    if line[3].startswith(" Defective"):
                        reactomid = line[1]
                        print(reactomid)
                        # print(line)
                        break
                    elif namegene == line[3]:
                        reactomid = line[1]
                        print(reactomid)
                        # print(line)
                        break
                    elif namegene in line[5]:
                        reactomid = line[1]
                        print(reactomid)
                        # print(line)
                        break

                elif ensmbleid not in line:
                    reactomid=''


            url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id='
            out = open('out.txt', 'w')
            response = requests.get(url + ncbiid)
            out.write(response.text)
            out.close()
            out = open('out.txt', 'r').readlines()
            for line in out:
                print(line)
            for line in out:
                if '<Name>' in line:
                    name = line.strip().replace('<Name>', '').replace('</Name>', '')
                if '<NomenclatureName>' in line:
                    nomenclature = line.strip().replace('<NomenclatureName>', '').replace('</NomenclatureName>', '')
                if '<Description>' in line:
                    description = line.strip().replace('<Description>', '').replace('</Description>', '')
                if '<MapLocation>' in line:
                    maplocation = line.strip().replace('<MapLocation>', '').replace('</MapLocation>', '')
                if '<OtherAliases>' in line:
                    otherAliases = line.strip().replace('<OtherAliases>', '').replace('</OtherAliases>', '')
                if '<ExonCount>' in line:
                    exoncount = line.strip().replace('<ExonCount>', '').replace('</ExonCount>', '')
                if '<GeneWeight>' in line:
                    geneweight = line.strip().replace('<GeneWeight>', '').replace('</GeneWeight>', '')
                if '<Summary>' in line:
                    summary = line.strip().replace('<Summary>', '').replace('</Summary>', '')
                if '<ScientificName>' in line:
                    sciname = line.strip().replace('<ScientificName>', '').replace('</ScientificName>', '')
                if '<OtherDesignations>' in line:
                    designation = line.strip().replace('<OtherDesignations>', '').replace('</OtherDesignations>', '')
                if '<ChrLoc>' in line:
                    chrloc = line.strip().replace('<ChrLoc>', '').replace('</ChrLoc>', '')
                if '<ChrAccVer>' in line:
                    chracc = line.strip().replace('<ChrAccVer>', '').replace('</ChrAccVer>', '')
                if '<ChrStart>' in line:
                    chrstart = line.strip().replace('<ChrStart>', '').replace('</ChrStart>', '')
                if '<ChrStop>' in line:
                    chrend = line.strip().replace('<ChrStop>', '').replace('</ChrStop>', '')

    genes = zip(PharmGKB_Accession_Id, NCBI_Gene_ID, HGNC_ID, Ensembl_Id, Name, Symbol, Alternate_Names,
                Alternate_Symbols, Is_VIP, Has_Variant_Annotation, Cross_references, Has_CPIC_Dosing_Guideline,
                Chromosome, Chromosomal_Start_GRCh37_p13, Chromosomal_Stop_GRCh37_p13, Chromosomal_Start_GRCh38_p7,
                Chromosomal_Stop_GRCh38_p7)

    print(summary,'********************')
    print(ncbiid)

    return render(request, 'genes.html',context={'genes': genes, 'apiname': name, 'description': description, 'maplocation': maplocation,
                                                    'otherAliases': otherAliases,'summary':summary, 'exoncount': exoncount,
                                                    'geneweight': geneweight, 'sciname': sciname, 'designation': designation,
                                                    'chrloc': chrloc, 'chracc': chracc, 'chrstart': chrstart, 'chrend': chrend,
                                                    'nomenclature': nomenclature,
                                                    'ensid': ensmbleid, 'ncbiid': ncbiid, 'reactomid': reactomid})




def table(request):

    PharmGKB_Accession_Id = []
    Name = []
    Generic_Names = []
    Trade_Names = []
    Brand_Mixtures = []
    Type = []
    Cross_references = []
    SMILES = []
    InChI = []
    Dosing_Guideline = []
    External_Vocabulary = []
    Clinical_Annotation_Count = []
    Variant_Annotation_Count = []
    Pathway_Count = []
    VIP_Count = []
    Dosing_Guideline_Sources = []
    Top_Clinical_Annotation_Level = []
    Top_FDA_Label_Testing_Level = []
    Top_Any_Drug_Label_Testing_Level = []
    Label_Has_Dosing_Info = []
    Has_Rx_Annotation = []
    num=[]
    inpo =  request.GET.get('drug' ,'M')

    file = open('staticfiles/chemicals.tsv',
                'r').readlines()
    for line in file:
        line = re.split('\t', line.strip())
        if inpo in line[1]:
            try:
                PharmGKB_Accession_Id.append(line[0])
                Name.append(line[1])
                Generic_Names.append(line[2])
                Trade_Names.append(line[3])
                Brand_Mixtures.append(line[4])
                Type.append(line[5])
                Cross_references.append(line[6])
                SMILES.append(line[7])
                InChI.append(line[8])
                Dosing_Guideline.append(line[9])
                External_Vocabulary.append(line[10])
                Clinical_Annotation_Count.append(line[11])
                Variant_Annotation_Count.append(line[12])
                Pathway_Count.append(line[13])
                VIP_Count.append(line[14])
                Dosing_Guideline_Sources.append(line[15])
                Top_Clinical_Annotation_Level.append(line[16])
                Top_FDA_Label_Testing_Level.append(line[17])
                Top_Any_Drug_Label_Testing_Level.append(line[18])
                Label_Has_Dosing_Info.append(line[19])
                Has_Rx_Annotation.append(line[20])
            except:
                continue

    for val in range(100):
        num.append(str(val))

    drug = zip(PharmGKB_Accession_Id, Name, Generic_Names, Trade_Names, Brand_Mixtures, Type, Cross_references, SMILES,
               InChI, Dosing_Guideline, External_Vocabulary, Clinical_Annotation_Count, Variant_Annotation_Count,
               Pathway_Count, VIP_Count, Dosing_Guideline_Sources, Top_Clinical_Annotation_Level,
               Top_FDA_Label_Testing_Level, Top_Any_Drug_Label_Testing_Level, Label_Has_Dosing_Info, Has_Rx_Annotation , num)
    print(Name)
    return render(request, 'drug.html', context={'drug':drug})


def phylo(request):
    snpo = request.GET.get('phylo', None)
    if snpo == None:
        return render(request, 'phy.html')
    else:
        selection = request.GET.get('selection')
        if selection == None:
            raise Exception('You Did not Choose input type')
        elif selection == '1':
            snp = []
            if True:
                Entrez.email = "A.N.Other@example.com"  # Always tell NCBI who you are
                snp_seq = open("staticfiles/Motif_in.fasta", 'w')
                snpo = re.split('\s', snpo)
                for val in snpo:
                    snp.append(val)
                    print(val)
                print(snp)

                seqlenghtes = []
                idslengthes = []
                httpseq=[]
                httpid=[]
                for id in snp:
                    entrez_id = id.strip()
                    handle = Entrez.efetch(db="snp", id=entrez_id, rettype="fasta", retmode="fasta")
                    record = SeqIO.parse(handle, "fasta")
                    for rec in record:
                        print(rec)
                        httpseq.append(rec.seq[0:20]+'...')
                        httpid.append(rec.id.replace('gnl|dbSNP|',''))
                        idslengthes.append(len(rec.id))
                        seqlenghtes.append(len(rec.seq))
                        snp_seq.write('>%s\n%s\n\n' % (rec.id, rec.seq))

                print(idslengthes)
                print(seqlenghtes)
                snp_seq.close()
                x = max(int(s) for s in idslengthes)
                seqlen = min(int(s) for s in seqlenghtes)
                from Bio import SeqIO
                records = SeqIO.parse("staticfiles/Motif_in.fasta", "fasta")
                lens = []
                lens2 = []
                file = open("phylo.phy", 'w')
                file.write('       %s      %s\n' % (len(idslengthes), int(seqlen)))

                for i, record in enumerate(records):
                    ids = record.id
                    sequence = record.seq[0:int(seqlen)]
                    lens.append(record.id)
                    lens2.append(record.seq)
                    line = "%s   %s" % (ids, sequence)

                    print("")

                    lengthmax = len(max(lens, key=len))
                    lengthmin = len(min(lens, key=len))

                    idlen = len(lens) - 1

                    if idlen < int(x):
                        # print(len(record.id), x)
                        add = int(x) - len(record.id)
                        if i in range(0, 10):
                            iden = str(i) + '-' + record.id + (add * '-') + '--'
                        elif i in range(10, 100):
                            iden = str(i) + '-' + record.id + (add * '-') + '-'
                        else:
                            iden = str(i) + '-' + record.id + (add * '-')
                        print(iden, sequence)
                        iden = iden

                        file.write('%s  %s\n' % (iden.replace('gnl|dbSNP|', ''), sequence))
                    else:

                        # print(len(record.id), x)
                        add = int(x) - len(record.id)
                        if i in range(0, 10):
                            iden = str(i) + '-' + record.id + (add * '-') + '--'
                        elif i in range(10, 100):
                            iden = str(i) + '-' + record.id + (add * '-') + '-'
                        else:
                            iden = str(i) + '-' + record.id + (add * '-')
                        print(iden, sequence)
                        file.write('%s    %s\n' % (iden.replace('gnl|dbSNP|', ''), sequence))
                file.close()
                alignment = []
                # Read the sequences and align
                aln = AlignIO.read('phylo.phy', 'phylip')
                for align in aln:
                    alignment.append(align)

                # Print the alignment
                print(aln)

                matrix = []
                # Calculate the distance matrix
                calculator = DistanceCalculator('identity')
                dm = calculator.get_distance(aln)

                # Print the distance Matrix
                print('\nDistance Matrix\n===================')
                print(dm)
                for val in dm:
                    matrix.append(val)

                # Construct the phylogenetic tree using UPGMA algorithm
                constructor = DistanceTreeConstructor()
                tree = constructor.upgma(dm)
                



                from Bio.Seq import Seq
                from Bio import motifs
                import os
                list2 = []
                file = open('staticfiles/Motif_in.fasta', 'r').readlines()
                out = open('staticfiles/Motif_out.fasta', 'w')
                for line in file:
                    line = line.replace('N', '').replace('D','').replace('W','').replace('K','').replace('t','T').replace('a','A').replace('c','C').replace('g','G').replace('R','')
                    out.write(line)
                out.close()

                for record in SeqIO.parse('staticfiles/Motif_out.fasta',
                                          "fasta"):
                    print(record.seq)
                    list2.append(record.seq[0:400])

                m = motifs.create(list2)
                m.weblogo("staticfiles/MOTIF2.svg", format='SVG', stacks_per_line=100, stack_width='large',
                          color_scheme='color_auto', show_fineprint=True, scale_width=True)
                consensus=m.consensus



                # Print the phylogenetic tree in the terminal
                print('\nPhylogenetic Tree\n===================')
                Phylo.draw_ascii(tree)
                # Draw the phylogenetic tree
                Phylo.draw(tree)

                # Print the phylogenetic tree in the terminal
                print('\nPhylogenetic Tree\n===================')
                Phylo.draw_ascii(tree)

                # Draw the phylogenetic tree
                def get_label(leaf):
                    return leaf.name

                tree.ladderize()
                Phylo.draw(tree, label_func=get_label, do_show=False)
                pylab.axis('on')
                pylab.savefig('staticfiles/tree2.svg', format='svg', bbox_inches='tight', dpi=300)

                all_data=zip(httpid,httpseq)
    return render(request, 'phy.html', context={'inpo': alignment,'all_data':all_data,'cons':consensus})


def phylogene(request):

    if True:
        snp=[]
        Entrez.email ="A.N.Other@example.com"  # Always tell NCBI who you are
        import re
        choice = request.GET.get('phylo', 'C9orf72')
        print(choice)
        genesfile = open(
            'staticfiles/var_pheno_ann.tsv',
            'r').readlines()
        for line in genesfile:
            line = re.split('\t', line)
            if choice in line[2]:
                print(line[1])
                snp.append(line[1])

        snp_seq = open("SNP_seq.txt", 'w')
        # snp = re.split('\s', str(snp))
        print(snp)

        lenghtes = []
        idslengthes = []

        for id in snp:
            entrez_id = id.strip()
            handle = Entrez.efetch(db="snp", id=entrez_id, rettype="fasta", retmode="fasta")
            record = SeqIO.parse(handle, "fasta")
            for rec in record:
                SeqIO.write(rec, "" + entrez_id + ".fasta", "fasta")
            # print(record.id, record.seq)
            snp_seq.write('>%s\n%s\n' % (rec.id, rec.seq))
            idslengthes.append(len(rec.id))

        print(idslengthes)
        print(len(idslengthes))
        snp_seq.close()
        x = max(int(s) for s in idslengthes)
        records = SeqIO.parse("SNP_seq.txt", "fasta")
        lens = []
        lens2 = []
        file = open("phylo.phy", 'w')
        file.write('       %s      %s\n'%(len(idslengthes),'100'))


        for i,record in enumerate(records):
            ids = record.id
            sequence = record.seq[0:100]
            lens.append(record.id)
            lens2.append(record.seq)
            line = "%s   %s" % (ids, sequence)


            print("")

            lengthmax = len(max(lens, key=len))
            lengthmin = len(min(lens, key=len))

            idlen = len(lens) - 1

            if idlen < int(x):
                # print(len(record.id), x)
                add = int(x) - len(record.id)
                if i in range(0, 10):
                    iden = str(i) + '-' + record.id + (add * '-') + '-'
                else:
                    iden = str(i) + '-' + record.id + (add * '-')
                print(iden, sequence[0:100])
                file.write('%s    %s\n' % (iden, sequence[0:100]))
            else:

                # print(len(record.id), x)
                add = int(x) - len(record.id)
                if i in range(0, 10):
                    iden = str(i) + '-' + record.id + (add * '-') + '-'
                else:
                    iden = str(i) + '-' + record.id + (add * '-')
                print(iden, sequence[0:100])
                file.write('%s    %s\n' % (iden, sequence[0:100]))
        file.close()

        if choice == '':
            return render(request, 'phygene.html')
        else:

            alignment = []
            # Read the sequences and align
            aln = AlignIO.read('phylo.phy', 'phylip')
            for align in aln:
                alignment.append(align)

            # Print the alignment
            print(aln)

            matrix = []
            # Calculate the distance matrix
            calculator = DistanceCalculator('identity')
            dm = calculator.get_distance(aln)

            # Print the distance Matrix
            print('\nDistance Matrix\n===================')
            print(dm)
            for val in dm:
                matrix.append(val)

            # Construct the phylogenetic tree using UPGMA algorithm
            constructor = DistanceTreeConstructor()
            tree = constructor.upgma(dm)

            # Print the phylogenetic tree in the terminal
            print('\nPhylogenetic Tree\n===================')
            Phylo.draw_ascii(tree)

            # Draw the phylogenetic tree
            Phylo.draw(tree)
            return render(request, 'phygene.html')

def pathway(request):
    gene1 = request.GET.get('gene1','CFTR')
    gene2 = request.GET.get('gene2','VCP')
    organism = request.GET.get('organism','Mus musculus')
    import re
    inpo = []
    inpo.append(gene1)
    inpo.append(gene2)
    inpo.append(organism)
    print(inpo)
    file = open('staticfiles/Pathway_NCBI.tsv',
                'r').readlines()
    for line in file:
        line = re.split('\t', line)
        if inpo[0] in line[3]:
            try:
                if inpo[1] in line[3]:
                    try:
                        if inpo[2] in line[5]:
                            print('')
                            reid=line[1]
                    except IndexError:
                        print('')
                        reid = line[1]
            except IndexError:
                print('')
                reid = line[1]

    return render(request, 'pathway.html', context={'id':reid})

def disease_database(request):
    import re
    disease_name_list = []
    snp_id_list = []
    chemical_name_list = []
    clinical_annotation_list = []
    gene_name_list=[]
    dis = request.GET.get('disease', 'Diabetes')
    disease_file = open(
        'staticfiles/clinical_ann_metadata.tsv',
        'r').readlines()
    for disease in disease_file:
        disease = re.split('\t', disease)
        if dis in disease[12]:
            disease[12] = re.split('\s', disease[12])
            disease[11] = re.split('\s', disease[11])
            disease[2]= re.split('\s' , disease[2])
            snp_id = disease[1]
            disease_name = disease[12][0].replace('"','')
            chemical_name = disease[11][0]
            clinical_annotation = disease[8]
            gene_name=disease[2][0]
            gene_name_list.append(gene_name)
            clinical_annotation_list.append(clinical_annotation)
            snp_id_list.append(snp_id)
            chemical_name_list.append(chemical_name)
            disease_name_list.append(disease_name)
            print(snp_id, " ", chemical_name, " ", disease_name, " ", clinical_annotation)
    all_data=zip(disease_name_list,clinical_annotation_list,snp_id_list,chemical_name_list,gene_name_list)
    return render(request, 'disease_database.html' , context={'all_data':all_data})

def disease_catcher(request):
    import re
    chromosome = '17'
    gene = 'BRCA1'
    q = 'AGCACTGACGCCAGATGTAATTTTGGTTTGAAGTTCTCAAGACTGGCAAGGAGCTTCACCACAGTTTAAGCAAGCTCAAGCCTAAGGAAGATAGGATGTGAAAAACTGGGTGGTCACATGGGACTCTCAGAAGTCTCAATGCATTTCGACTCCTTTCATTATTTTAAAAGAGACCAATCTGCCTACAAAAGGGAAAACAAAAGAAAAGCAGAAGTTACTTCACTTTAAGGAACTTCCATAAGTGATTGTCAGAGCATTTAATTTTTCTCTGGAAAAGTTGCTGAGCTGTCATTCAGTTTAAGGGAAGGAAGTTGGATAGGGCAGAAAAGAGATTCTTTGAATCAAGGAGAGAGGAGGCCAGACTGGTCACAGGTGTGGTTAAGACCCCTATCTAAAGAAGCCCTGGGAGAAAGCAGACACCTCACTGCTTTCTGAGGCGGGGTAAGAAAGCCTGAGGAGGGTCATAGTGATACCAGGGTCTGTCTTTCACTAGGAAAACTCAAGGCTAAATCTGGACTCTGTGTGACAGAACGAGGCCGGCAAGTGAGCCACACAGATCTTTCTGGGCTAAAATATAGGGATGCTCACCTTGTTCAGGTCTGGGATCCCAGGCTCCTCAGCCGCACGAGGGCGAGACTTCTGGAGTGCGTCTTTGGGGAGAAGGGCCTCTAAGGGCAGAGTGGGCTGGGCTGGGCTGAGCTGAGCTGGGCGGTGGTGGGGAGAAGCACCGCCTGGGGCCGGCAACCGGCCGGAGTCCTAGGGGAACGAATGTGCCGCACTGGAGAGTGGGGACAGTAAGGACGCTGACCTCCCCGAAGGGCCAAGCCGATTCCTGGGCTACATGGGGGCTCTGCCATGGGCCGACCGGGTGGAAGCGGACGGAGAAGGTGTGACTTCATTAAGGGAGGATTTCAAGTCCCACCACACACCCCTAGAGGCCGTTAGTCAAGCTTGGTACCCGTGCGGGGGGCGCCTGTGTGTTGGGCGAGGAGGCGGCAACAGGAGGTGTCTTGACCTCCGAAAGGTAAAGTAGGTATAAAGGGTCTAGGGAGGGCCTTACCCGCGGACGACTGCAGCGACAGAAACGGAAGAGAGTGGAATCGAAAATGCGGGAGCGAACACTCCCCTCAAGCTTTGACTCAGTGCCCAGGGAAGTCCCGCCCCAGCCCCTAAGGGCCCACCCCCCCGTCCCCCGGTCCGCCCCTTCGCCACGCAGCAGCCTGAGACGCTGCCAGCTCTTCAAACGTTACCCGAGGGCCCGGACGCTGCTTCCGCAAGGACTAGCCAGTCCGGGTAGGAATTTCCAAGAACATTTATTTCACCATTCTTTAGGTGGTTATGTACAAGTACATACATAACCAGGAGAGCAACCCGGGCCCTCCTGGCCCTCGATTATAAAACTCCGATGGTCAAGACAGGTCCCACCCATCAACCCAATTACACCCTGCCGCTCCACCTTCTACCCGGAAGCCCCGATTAAAACGTTTTAAATAATGAACGGCCTGGCTGGCTTGCAGTTTGAATCCTGGTTTATTTAAGTAGCATTCAGGAGTTATACTGTAAATAAACATCATTATGGCTAATGACTAATTAAAATCACCATTCATTACAAATACCAACCAGGTTAGATAATACTATTAACTATCTTTTCCCCTTCGGTCTGGGGTTCATATTCTTCCTTTTCTGGAAGTCATTTTCATGGATAAAATATTTAGCATTTTTTGTTACCTCTTCCTTCACCAACTGTACTCAACACACTTCTGTAAACAGAGTGAGAACCTGGGTTCATGGACAAGCTCTTTTCCACTTATCTTCGGGTTAAACCAAAACTTTTTCAGCAACTTTGCCCCTGTCCAAGTTTTGCAGAACACCAGCTGCCTCTTCCTGTTCTAGGCATTCTAGGGCTAGATACCTGTGAACGGACTTCAAGCATTAAAGAAGTCATCAGCATTTCCTATACCTCACACACACTCATTTGCCCACCCTGAACTTCACTGACATTAATTTAGGAGCCTGCCTTTTGCCTTTCACATGTTTCTATTTCCCCCACTGTACTGTAAAATAAGATCTGCCCTGACATTAAGAAAAAAGGCAACTGTTAAATTTACTTTTCCCTTTAAATAACGGGATTTTTTCCTTTTCTATAAGATTAGCAAAATATTTCACAATACCTATATTTTTTCCCTTAAATAAAAATGTATGCAAAATTACTATAATAGACATTGCTCAATAGGTGCTGGTAGGTTTATCAGCACACTCTTCCTCCTTTCATTGAAATGCGGAGATTCCCTTGAAGTTGGATACTCTTGAATGCACCATGGAGTCTTTAGAAAGCAACGGTGCCAATTATTAAGTATTGCAGGAAAATTACTTTCTCCTTGAGAGAGGAAGAGTAACAACATTTCTATCATGCTCTGGGTGTTATTTCAATGCCAAATAACTATGAAGCTTGAGGTGGATTACAGAAATTATGAGAATAAGAATTATGAGAATAATTCTTAAAAACAACTTCCGAGGTTAATGTCTTTCAAAACCTGAAAATTAACAGATATCCCCATAAAAAATCCACATTAGCCAAAGAAAAATACTCAGTTTCAATAAAAGCAGAAAGAAAAATGGTTCTATATGAAGTTGAAGTTACACTTGTCAATGCACTAAAATAGTTTTATCCAAATCAACTGAAAAACTGTCTTCATACTATTCATCCAATGTCTAAAGTGAAAATTTAAGTTCTGGTCTTTAGTTTGCTCCCCAAAATACTGTAGCAGGTATTGTCTTGGTAACTCAACGATGGGCTCTATACAGATTCATCAGATTAAAAATAAGCAAGGGCTTCTATCCCATACCCCAATGACAGAATAAAGATCATCTCTGAGCAGCAGGCAGTTATTTAATACAAGGTCACTCCTCAATAGCGTTGGCTGTACCAGTCGTTGTTTTAACTGAAGAAGTTCTGTCACAACCTGCAGGATATTGTAATTGTGTTTCTTCAGCAGCCGTAGGTTCAGCTGCCTGTCACAGAATCCCATTTCAAAGAGATGGGCCATCAGGGCTGCTGTCTGATCTTCAGAAATTATTGGCTGTGTAGAGACAAGCAAAGCAAAAAAAAAAAAAAAAATTCTGTATAATTTAAATGTTTATTTTATACATTTGTATCATCTTGCAAGCATACTAGTCTTTCCAGATTTATTCAAAGTACTAGGTAATTAGATTCAAGGCTATGGCACTCCTGTGTACTTGATACCAACGGTTCGGACTGGGAAGGCCCTATTTTAGATTAGTTATCTGTGTACTTGAAGGAGGAGAAAAGGGAAACCGGGGAGCTAATGAACTGTACTTAGAACCAGGAGAATGCTTCACTCGGGGCTTTCAGGAAGAAGCATTTTTTGGCCAGATGAAGGAATGACAGCAAGCCCCAGAGACCCTGTGAACATAAGTGTAATTGTGATGATGGGGGAAGTATGAGCAAGATAAAAGTGTTACTGAAAGATGCTATGGGAAATTTTGTCAGATATGCAGAGGAGCTGGACAAAAAAAACAGCTTGGGAAATTACAACCTGGGAAGGCTTCTTGATTTTTTCTCCTCCAAAAGAATAGGTAAGAGGGGCTGGGTGCGGTGGCTCACACCTGTAATCCCAGCACTTTGGGAGGCTGAGGCGGGCGGATCATGAGGTCAGGAGATTGAGACCATCCTGGCTAACACAGTGAAACCCCATCTCTACTAAAAATACAAAAAATTAGCCAGGCGTGGTGGCAGGCGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGCAGGACAACGGCGTGAACCCAGGAGGTGGAGCTTGCAGTGAGCCGAGATCGCGCCATTGCACTCCAGCCTGGGTGACAGAGCGAGACTCCATCTCAAAAAAAAAAAAACAAGAATAGGTAAGAGGAAGTTAACATTCATTCATTCAAGCCCTGTGGTGGAAGCATTACCTTCCTAAATCCTTAAACGCTTGTTTTGTTGTTGTTGTTGTTTTGTTTGTTTTTTGAGACAGGGTGTCACTCTGTTCCCCAGGCTGGTGTGCAGTGGTGTGAACACAGCTCATTGCAGCCTCAACCTCCTGGGCTCAAGTGATCCTCCTGCCTCAGCCTTCCAAATAGCTGGGACCACTGGTCTGCGCTACCACGCCCGGCAAATTTTTAAATTTTTTACTAGAGATGAGGTTTCCCTATGTTGCCCAGGCTGGTCTTGAACTTCTGGGCCCAAGTGATCCTTGAACTTTGGTCTCCCAAAATGCTGGGATTATAGATGTGAGCCACCTCGCCTGGCCTCTGGTTCTGTTATTATCCCAATTTCACAGACTGGGGATACTGAAACTGTGCAGCAGAAAGATTATTAACTTGGGAGGCAGATAGGCTTAGACTCAAACCCTAATCTTCCATTTACTTACCAACTGTGCTATCAATTCCTTTATGAACTGAAGATATATCATACCTCCTACTTTGCAGGGTTGTTGAAAAAGTATGTAAAGATGAGAACAGTATGCAAACATGTTAAGCATTTCCCTCCAAACCTGATTCCACTGTTCCATTCTACCTTCTTATTTCAGAATAGATTATCCATCCTCCTTTGCAGGTTAATGGAATGTACCTTGCATTGTATCTTATTACCTCTCACCACTTCAGAAATCTTTCTCCATCTCCAATCTGTCTTCTCATCTTTATCATTAACCTCTTTCACTGCTGGCTTGCCCATCTCTAAAGTCTCACCCAACTTTAAAAATGCTTCCTCTTGTCTTATGTATCCTAGCTTCTGCCTGTTGTAGCCTGTCTATTTGAAAGTGTGGCATACTTCTACTCTCAAATTCCTCACCACCAGTTCATATCTCTCTGTAGCATCTAACACCGTTGACCATGTCCCCCTTATTAAAACACCCTTCCTTAAGTTCCAGCATAGAACTCTCTCTCTTCATTTTCCTCCTAGTTCCCAGCTTTCTTTCAGAGTATTCCTTTCCCTGATCCTTAAAGGTTGC'
    blast_out = open('/home/peter/Desktop/BLAST-API-Implementation-master/examples/blast_out.txt', 'r').readlines()
    for result in blast_out:
        result = re.split('\t', result)
        try:
            qstrat = result[6]
            qend = result[7]
            sstart = result[8]
            send = result[9]
            # print(result,qstrat,qend,sstart,send)
            break
        except:
            pass

    lenghtes = open('staticfiles/Gene_Lengthes.txt',
                    'r').readlines()
    for line in lenghtes:
        line = re.split('\t', line)
        if line[7] == gene:
            if line[2] == chromosome:
                gene_start = int(line[3]) + int(qstrat) - 1
                gene_end = int(line[3]) + int(qend)
                # print(line)
                # print('gene start from',gene_start,'and end on',gene_end)
                break
            break

    # lengthes_2=open('/home/peter/Desktop/Biotechnology/Pharmosome/Pharmosome/PharmGBK/genes/genes.tsv',)

    snp_count = set()
    before_snp = set()
    after_snp = set()
    snp_list = set()
    rs_list = set()
    snp_file = open(
        '/home/peter/Desktop/Biotechnology/Pharmosome/Pharmosome/PharmGBK/SNP_LOC_POS/snp_out_BRCA1_CHR(1).txt',
        'r').readlines()
    for snp in snp_file:
        snp = re.split('\s', snp)
        # print(snp)
        snp_name = snp[0].replace('"name":"', '').replace('"', '')
        snp_pos = snp[2].replace('"pos":', '')
        snp_ref = snp[4].replace('"ref":"', '').replace('"', '')
        snp_alt = snp[6].replace('"alt":"', '').replace('"', '')
        snp_chr = snp[8].replace('"chr":"', '').replace('"', '')
        # print(snp_name,snp_pos,snp_ref,snp_alt,snp_chr)
        if gene_start <= int(snp_pos) <= gene_end:
            # print(int(snp_pos)-int(gene_start))
            query_snp = q[int(snp_pos) - int(gene_start):int(snp_pos) - int(gene_start) + 1]
            if query_snp == snp_ref:
                pass
            elif query_snp == snp_alt:
                snp_info = 'SNP', snp_name, 'found in position', int(snp_pos) - int(gene_start), snp_ref, snp_alt
                snp_count.add(snp_name)
                before_snp.add(q[:int(snp_pos)])
                snp_list.add('[%s>%s/%s]' % (snp_ref, snp_alt, snp_name))
                rs_list.add(snp_name)
                after_snp.add(q[int(snp_pos):])
    all_snp_data = zip(before_snp, snp_list, rs_list, after_snp)

    print('found', len(snp_count), 'SNP which are the following:')

    for s in snp_count:
        end = int(s[3])
        print(s)
        q = q[:end] + '[%s>%s/%s]' % (s[4], s[5], s[1]) + q[end:]
    #     before_snp.append(q[:end])
    #     snp_list.append('[%s>%s/%s]' % (s[4], s[5], s[1]))
    #     after_snp.append(q[end:])
    #     rs_list.append(s[1])
    # print(q)
    # print(before_snp)
    # print(snp_list)
    # print(rs_list)
    # print(after_snp)
    # for before_snp, snp, rs, after_snp in all_snp_data:
    #     print('')
    #     print(before_snp,snp,rs,after_snp)

    return render(request, 'disease_catcher.html',context={'q':q,'snp_count':snp_count,'snplen':len(snp_count),'rs':all_snp_data})





def snp_range(request):
    dis = request.GET.get('disease', None)
    selection = request.GET.get('selection')
    all_snp_list=[]
    if selection=='1':

        print(selection)
        from_gene=open('staticfiles/var_pheno_ann.tsv','r').readlines()
        for line in from_gene:
            line=re.split('\t',line)
            if dis in line[2]:
                print(line[1],'**********')
                all_snp_list.append(line[1])

    elif selection == '2':
        from_gene = open(
            'staticfiles/var_pheno_ann.tsv',
            'r').readlines()
        for line in from_gene:
            line = re.split('\t', line)
            if dis in line[11]:
                print(line[1])
                all_snp_list.append(line[1])

    elif selection == '3':
        from_gene = open(
            'staticfiles/var_pheno_ann.tsv',
            'r').readlines()
        for line in from_gene:
            line = re.split('\t', line)
            if dis in line[3]:
                print(line[1])
                all_snp_list.append(line[1])


    elif selection == '4':
        from_gene = open(
            'staticfiles/clinical_ann_metadata.tsv',
            'r').readlines()
        for line in from_gene:
            line = re.split('\t', line)
            if dis in line[12]:
                print(line[1])
                all_snp_list.append(line[1])


    print(selection,'-------')

    print(dis)
    print(all_snp_list)
    length=len(all_snp_list)


    return render(request, 'SNP_range.html', context={'all_snp_list':all_snp_list,'length':length,'dis':dis})






