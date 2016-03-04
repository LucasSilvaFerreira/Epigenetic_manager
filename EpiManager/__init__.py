__author__ = 'lucas'
# -*- coding: utf-8 -*-
import sys
import pybedtools
import re
from subprocess import Popen, PIPE
import os
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd


def open_file_fast(nome_do_arquivo, delimitador='\n', tabulador_interno='none'):
    '''Abre um arquivo e retorna um array contendo cada linha lida como um index desse array.
    Args:
    nome_nome_do_arquivo (string)= caminho do arquivo a ser aberto.
    delimitador (string) = delimitador que define o final de cada linha dentro do arquivo. defaut ('\n')
    tabulador_interno (string) =  Usado para utilizar um delimitador a nivel de linha. Criando um array com a linha splitada
    pelo delimitador escolhido. defaut('none')
                                                t (tab)
                                                s(espaço)
                                                ;(ponto e virgula)
                                                :(dois pontos)
                                                n (enter)
                                                '''

    if tabulador_interno == 'none':
        return [linha for linha in open(nome_do_arquivo, 'r').read().split(delimitador) if linha]
    elif tabulador_interno == 't':
        return [linha.split('\t') for linha in open(nome_do_arquivo, 'r').read().split(delimitador) if linha]
    elif tabulador_interno == 's':
        return [linha.split(' ') for linha in open(nome_do_arquivo, 'r').read().split(delimitador) if linha]
    elif tabulador_interno == ';':
        return [linha.split(';') for linha in open(nome_do_arquivo, 'r').read().split(delimitador)if linha]
    elif tabulador_interno == ':':
        return [linha.split(':') for linha in open(nome_do_arquivo, 'r').read().split(delimitador) if linha]
    elif tabulador_interno == 'n':
        return [linha.split('\n') for linha in open(nome_do_arquivo, 'r').read().split(delimitador) if linha]
    else:
        sys.stderr.write('Erro:\nEntre com um delimitador de linhas aceitavel: \n      t (tab) \n      s(espaço) \n'
                         '      ;(ponto e virgula) \n      :(dois pontos) \n      n(enter)\n')
        sys.exit(1)

def bed_flank(bed_string, down= 200, up= 200):
    '''Given a bed_string, return a string with the two flank regions
    Parameters:
            bed_string : string bed separated por tab
                String in bed format

            down/up : int defaut(1000)
                The count of nucleotides extended in boundaries to create the flanks


    returns: str
        string with the two flanks separated for '\n'

    '''

    br = bed_string.split('\t')

    br[2] =  str(int(br[1])-1)
    br[1] = str(int(br[1])-down)
    br='\t'.join(br)

    bl = bed_string.split('\t')
    bl[1] =  str(int(bl[2])+1)
    bl[2] = str(int(bl[1])+up)
    bl='\t'.join(bl)

    return '\n'.join([br, bl])

def bed_start_site_to_tss(bed_string, down = 1000, up = 1000, not_reverse= True):
    '''dado uma coordenada *bed_string* retorna aquela string com o inicio e fim da coordenada sendo a extensão
        do start site dessa sequencia dada
        Parametros:
        -------------------
        bed_string : string bed separada por tab
            String no formato bed
        down/up : int defaut(1000)
            Quantidade de nuceotideos que serão extendidas partindo do start site do bed dado
        not_reverse : true
            Deixa os valores com coordenadas menores na primeria posição mesmo que a fita seja negativa
        '''
    saida_str = ''
    #print bed_string

    print type(bed_string), bed_string


    if "\t-" in bed_string:
        bed_string_tabulada =  bed_string.split('\t')
        ref_coord = int(bed_string_tabulada[2])
        n_start= ref_coord + up
        n_end= ref_coord - down
        bed_string_tabulada[1]= str(n_start)
        bed_string_tabulada[2]= str(n_end)
        #print (n_start)-(n_end)
        saida_str = '\t'.join(bed_string_tabulada)
    else:
        bed_string_tabulada =  bed_string.split('\t')
        ref_coord = int(bed_string_tabulada[1])
        n_start= ref_coord - down
        n_end= ref_coord + up
        bed_string_tabulada[1]= str(n_start)
        bed_string_tabulada[2]= str(n_end)
        saida_str= '\t'.join(bed_string_tabulada)

    saida_splitada= saida_str.split('\t')

    if int(saida_splitada[1]) >= 0 and  int(saida_splitada[2]) >= 0 : #caso nao exista valores negativos
        if not_reverse:
            if int(saida_splitada[1]) > int(saida_splitada[2]):
                temp1= str(saida_splitada[1])
                temp2= str(saida_splitada[2])
                saida_splitada[1]=temp2
                saida_splitada[2]=temp1
                return '\t'.join(saida_splitada)
            else:
                return saida_str
        else: #caso não queira o dado modificado para reversed
            return saida_str
    else:
        sys.stderr.write('valor proximo da borda do cromossomo: '+ saida_splitada[1]+'\n')
        return None

def bed_get_midle_point(bed_string):
    '''Given a *bed_string* respecting the 5' -> 3' coordinates
    returns a relative midpoint.
    Parameters:
    -------------------
    bed_string : str
        bed_string separated  by tab
    -------------------

    returns list:

    return [0]:
    string with mid point and the coords with midpoint +1
    ex: mid_point 100
    return [1]:

        ex: mid_point 100
        chr 100 101

    '''
    bed_string = bed_string.split('\t')
    #print bed_string

    midpoint = int(bed_string[2])-((int(bed_string[2])-int(bed_string[1]))/2)
    print 'midpoint', midpoint, bed_string[2], bed_string[1]
    second_return = bed_string
    print 'saida antes de processar', second_return
    second_return[1] = str(midpoint)
    second_return[2] = str(midpoint+1)
    print 'saida', '\t'.join(second_return)
    return str(int(bed_string[1]) + midpoint), '\t'.join(second_return)

def plot_bed_to_closest_bed_mark(bed_gene_file, bed_histone_mark, color='g' , plot_name_save=None, label=None, range_plot=None, show=False):
    '''Given a *bed_gene_file* coord file  and a *bed_histone_mark*, generates a distance density between bedfile and mark midle
    point mark
    Returns
    --------
        A bedtool object containing the distance in last row between *bed_gene_file* and *bed_histone_*mark'''
    import seaborn as sns
    import matplotlib.pyplot as plt


    # print ['\t'.join(bed_start_site_to_tss(bed_line).split('\t')[0:3])
    #                                           for bed_line in open_file_fast(bed_gene_file)
    #                                           if bed_line and not re.search('^#', bed_line)]


    bed_gene_file_array = []
    for bed_line_clear in [bed_line for bed_line in open_file_fast(bed_gene_file) if bed_line and not re.search('^#', bed_line)]:
        if bed_start_site_to_tss(bed_line_clear) != None:
            #print bed_line
            bed_gene_file_array.append('\t'.join(bed_start_site_to_tss(bed_line_clear).split('\t')[0:3]))


    tss_point =  pybedtools.BedTool('\n'.join(bed_gene_file_array), from_string=True)


    mark_midle_point = pybedtools.BedTool('\n'.join([ '{}\t{}\t{}'.format(tss_midle.split('\t')[0], str(bed_get_midle_point(tss_midle)), str(int(bed_get_midle_point(tss_midle))+1))
                                                     for tss_midle in open_file_fast(bed_histone_mark) if tss_midle and not  re.search('^#', tss_midle)]), from_string=True)


    #print tss_point[0:3]
    #print mark_midle_point[0:3]
    closest_feature = tss_point.sort().closest(mark_midle_point.sort(), D='a')

    data_plot = [distance[6] for distance in closest_feature]

    sns.set(style="white", palette="muted")


    sns.distplot(data_plot, hist=False, color=color, kde_kws={"shade": False, 'lw': 3},  label=label  )
    if range_plot:
        plt.xlim(-range_plot, range_plot)
    if show:
        plt.show()
    if plot_name_save:
        plt.savefig(plot_name_save+'.eps', format='eps')
        plt.close()
    return closest_feature

def bam_read_count(bamfile):

    """ Return a tuple of the number of mapped and unmapped reads in a bam file """
    p = Popen(['bamutils', 'stats', bamfile], stdout=PIPE)
    mapped = 0
    unmapped = 0
    mapped = 0
    ref_line = False
    for line in p.stdout:
        if ref_line and len(line)>2:
            #print line
            mapped+= int(line.replace('\t\n', '').split('\t')[1])
        if 'Reference' in line:
            ref_line= True
    return mapped

def __normalization__(genomecov_file, factor):
    '''Given a factor and genomecov_file makes a normalization of bam depth'''
    return_normalization = []
    sys.stderr.write("Normalizing depth coverage..." + "\n")
    for x in open_file_fast(genomecov_file, tabulador_interno='t'):
        data = x
        #print x
        #print genomecov_file

        if float(data[-4]) != 0:
            data.append(str((float(data[-4])/int(factor))*1000000))
        else:
            data.append(str(0))
        return_normalization.append('\t'.join(data))
    return return_normalization

def bam_sort_and_index_alinged_number(bam_file):
    '''creates a sorted bam and keep original file
    creates a index file
    get a number of aligned reads
    return number of aligned reads; total reads, sorted file name'''
    #print 'teste', bam_file
    #print 'entrou funcao'
    out_name = 'sorted_'+bam_file.split('/')[-1]
    out_name_create = bam_file.split('/')
    out_name_create[-1] = out_name

    complete_out = "/".join(out_name_create)
    #print bam_file
    #print '<', complete_out, '>'
    comando_sorted = 'bamtools sort -mem 7000 -in  {bam_file} -out {bam_file_out}'.format(bam_file=bam_file, bam_file_out=complete_out)
    comando_index = 'bamtools index -in {sorted_in}'.format(sorted_in=complete_out)
    sys.stderr.write("---Processing---\n {}...\n".format(out_name))
    sys.stderr.write("Ordenando bam..." + "\n")
    #print comando_sorted
    os.system(comando_sorted)
    sys.stderr.write("criando index..." + "\n")
    os.system(comando_index)
    sys.stderr.write("Capturando informações de alinhamento do bam..." + "\n")
    aligned_size = bam_read_count(complete_out)
    #print 'aling infos:', aligned_size, type(aligned_size)

    return aligned_size, complete_out

def coverage_bins(bam_file, string_bed_bins, name, adjust_bed_in, to_string=False):
    '''Calculates the coverage given a bam_file and a bedfile
    Parameters:
         bam_file <str>: A bam file
         string_bed_bins <str> : A complete bed file inside a string
         name <str> : A string name for created file
         to_string <bool><false>: Returns a list with coverage
    '''
    temp_bed_name = '{}_bed_slice_temp.bed'.format(name)
    bed_temp = open(temp_bed_name, 'w')
    ### checking adjust_bed_bin parameters

    print adjust_bed_in
    if not any([x for x in adjust_bed_in.values()]):
        for x_line in [l for l in open_file_fast(string_bed_bins) if not re.search('^#', l)]:
            bed_temp.write(x_line + '\n')
        bed_temp.close()

    elif len([[f, x]for f, x in adjust_bed_in.iteritems() if x != False]) > 2:
        sys.stderr.write( 'Error so many args in adjustament parameter' + "\n")
        print [[f, x] for f, x in adjust_bed_in.iteritems() if x != False]
        sys.exit()

    elif adjust_bed_in['r'] != False:
        for x_line in [l for l in open_file_fast(string_bed_bins) if not re.search('^#', l)]:
            # print '-'*20
            # print x_line
            # print bed_get_midle_point(x_line)[1]
            # print bed_start_site_to_tss(bed_get_midle_point(x_line)[1], down=adjust_bed_in['r'][0], up=adjust_bed_in['r'][1])

            bed_temp.write(bed_start_site_to_tss(bed_get_midle_point(x_line)[1], down=adjust_bed_in['r'][0], up=adjust_bed_in['r'][1]) + '\n')
        bed_temp.close()


    elif adjust_bed_in['f'] != False:
        for x_line in [l for l in open_file_fast(string_bed_bins) if not re.search('^#', l)]:
            flanks =  bed_flank(bed_string=x_line,
                          down=adjust_bed_in['f'][0],
                          up=adjust_bed_in['f'][1])

            bed_temp.write(flanks+ '\n')
        bed_temp.close()


    #print 'bedtools coverage -abam {bam_file} -b temp/bed_slice_temp.bed'.format(bam_file=bam_file)
    os.system('bedtools coverage  -abam {bam_file} -b {bed_slice_temp}> {name_temp}_coverage_temp.bed'.format(bam_file=bam_file, name_temp=name, bed_slice_temp=temp_bed_name ))
    os.system('rm {}'.format(temp_bed_name))
    os.system('bedtools sort -i {name_temp}_coverage_temp.bed > {name}_sorted_coverage_temp.bed'.format(name=name, name_temp=name))
    os.system('rm {name_temp}_coverage_temp.bed'.format(name_temp=name))
    if to_string:
        os.system('rm {name}_+sorted_coverage_temp.bed'.format(name=name))
        return open_file_fast('{name}_sorted_coverage_temp.bed'.format(name=name))

    else:
        #os.system('rm {name}_+sorted_coverage_temp.bed'.format(name=name))
        return '{name}_sorted_coverage_temp.bed'.format(name=name)

def get_normalized_count_in_bed_region_given_a_mark(bam, bed_target, analysis_name, resize_bed=False, flank_bed=False):
    '''Given a bam file and a Bed_file computes the normalized (rpkm) inside bed regions
    Parameters:
        bam <str>: A bam file path
        bed_target <str>: A bed file path
        analysis_name <str> : A string with the unique mark_name
        resize_bed list<int> default<False>: The original bed will be resized with base in you center
                                ex: resize_bed = [100,500]
                                chr   100-start stop+500

        flank_bed list<int> default<False>: The original bed will be transformed in  2 flanks regions
                                ex: flank_bed = [100,500]
                                chr   100-start start-1
                                chr   stop+1    stop+500
    return:

        bed_target with the last column containing the normalized count value




    '''
    factor, sort_bam = bam_sort_and_index_alinged_number(bam)

    cov = coverage_bins(bam_file=sort_bam, string_bed_bins=bed_target, name=analysis_name, adjust_bed_in={'r':resize_bed, 'f':flank_bed})
    norma = __normalization__(cov, factor)
    os.system('rm {}'.format(cov))

    return norma

def get_several_normalized_count_in_bed_region_given_a_mark(hash_bam, bed_target,  r_bed=False, f_bed=False ):
    '''
    This is a batch version of  function  get_normalized_count_in_bed_region_given_a_mark()

    Parameters:
        hash_bam_bed <dict<str>> : Ex: {
                                        mark1_name_string:'path/to/file.mark1_bam',
                                        mark2_name_string:'path/to/file.mark2_bam',
                                        mark3_name_string:'path/to/file.mark3_bam'
                                        }
        bed_target  <str>: A bed_target file



        r_bed list<int> default<False>: The original bed will be resized with base in you center
                                ex: r_bed = [100,500]
                                chr   100-start stop+500

        f_bed list<int> default<False>: The original bed will be transformed in  2 flanks regions
                                ex: f_bed = [100,500]
                                chr   100-start start-1
                                chr   stop+1    stop+500





    '''

    all_results_dict = {name_analysis: get_normalized_count_in_bed_region_given_a_mark(bam=bam, bed_target=bed_target, analysis_name=name_analysis, resize_bed=r_bed, flank_bed=f_bed) for name_analysis, bam in hash_bam.iteritems()}



    df_to_box = pd.DataFrame({name:[float(value_count[-1]) for value_count in lista_out] for name, lista_out in all_results_dict.iteritems()})




    sns.boxplot(data=df_to_box, orient="h", palette="Set2")
    plt.xlim(-5,5)
    plt.show()

def get_groups_normalized_count_in_bed_region_given_a_mark(hash_in_groups_bam, bed_target_hash, r_bed=False, f_bed=False ):
    '''
    This is a group batch version of  function  get_normalized_count_in_bed_region_given_a_mark()

    Parameters:
        hash_in_groups_bam <dict<dict<str>>> : Ex: {

                                        'tss':{
                                        mark1_name_string:'path/to/file.mark1_bam',
                                        mark2_name_string:'path/to/file.mark2_bam',
                                        mark3_name_string:'path/to/file.mark3_bam'},


                                        'coding':{
                                        mark4_name_string:'path/to/file.mark5_bam',
                                        mark5_name_string:'path/to/file.mark6_bam',
                                        mark6_name_string:'path/to/file.mark7_bam'
                                        }





                                        }

        bed_target  dict<str>: Ex:{ 'tss':'path/to/tss.bed',
                                    'coding':'path/to/coding.bed'
                                 }


        r_bed list<int> default<False>: The original bed will be resized with base in you center
                                ex: r_bed = [100,500]
                                chr   100-start stop+500

        f_bed list<int> default<False>: The original bed will be transformed in  2 flanks regions
                                ex: f_bed = [100,500]
                                chr   100-start start-1
                                chr   stop+1    stop+500




    '''

    df_array_to_paste = []

    for bed_name, bed_file_group in  hash_in_groups_bam.iteritems():

        bed_target = bed_target_hash[bed_name]
        hash_group_bed = bed_file_group

        all_results_dict = {name_analysis: get_normalized_count_in_bed_region_given_a_mark(bam=bam, bed_target=bed_target, analysis_name=name_analysis, resize_bed=r_bed, flank_bed=f_bed) for name_analysis, bam in hash_group_bed.iteritems()}

        df_mix = []
        for name, lista_out in all_results_dict.iteritems():
            df_to_mix = pd.DataFrame({'Count':[float(value_count[-1]) for value_count in lista_out]})
            df_to_mix['Mark'] = name
            df_mix.append(df_to_mix)

        merged_marks_form_same_bed_source = pd.concat(df_mix)
        merged_marks_form_same_bed_source['Group'] = bed_name
        #[pd.DataFrame({'Count':[float(value_count[-1]) for value_count in lista_out]})for name, lista_out in all_results_dict.iteritems()]

        #df_to_box = pd.DataFrame()

        df_array_to_paste.append(merged_marks_form_same_bed_source)

    concact_df_array = pd.concat(df_array_to_paste)

    #print concact_df_array
    pd.DataFrame.to_csv(concact_df_array, 'tabela_grupo.txt' , sep='\t', index=False)
    sns.boxplot(x="Group", y="Count", hue="Mark", data=concact_df_array, palette="Set2")


    #sns.boxplot(data=df_to_box, orient="h", palette="Set2")
    plt.ylim(-5,5)
    plt.show()


def main():
    # plot_bed_to_closest_bed_mark(bed_gene_file='/work/users/lucassilva/chromhmm/sra_dir/peaks_hmm/AR_24/AR_24.bed',
    #                          bed_histone_mark='/work/users/lucassilva/chromhmm/sra_dir/peaks_hmm/PolII/PolII.bed',
    #                              show=True)
    #


    print len(get_normalized_count_in_bed_region_given_a_mark(bam= '/work/users/lucassilva/chromhmm/sra_dir/peaks_hmm/H3K4me1/H3K4me1.bam',
                                                    bed_target='/work/users/lucassilva/chromhmm/sra_dir/peaks_hmm/AR_24/AR_24.bed',
                                                    analysis_name='H3K4ME1', flank_bed=[10, 10]))

    print len(get_normalized_count_in_bed_region_given_a_mark(bam= '/work/users/lucassilva/chromhmm/sra_dir/peaks_hmm/H3K4me1/H3K4me1.bam',
                                                    bed_target='/work/users/lucassilva/chromhmm/sra_dir/peaks_hmm/AR_24/AR_24.bed',
                                                    analysis_name='H3K4ME1', resize_bed=[10, 10]))



    # get_several_normalized_count_in_bed_region_given_a_mark(hash_bam = {'H3K4ME1':'/work/users/lucassilva/chromhmm/sra_dir/peaks_hmm/H3K4me1/H3K4me1.bam',
    #                                                                     'H3K4MEteste':'/work/users/lucassilva/chromhmm/sra_dir/peaks_hmm/H3K4me1/H3K4me1.bam'
    #
    #                                                                     },
    #                                                 bed_target='/work/users/lucassilva/chromhmm/sra_dir/peaks_hmm/AR_24/AR_24.bed')
    #
    #


    bam_groups = {'exon':{'H3K4ME1':'/work/users/lucassilva/chromhmm/sra_dir/peaks_hmm/H3K4me1/H3K4me1.bam',
                              'H3K4ME3':'/work/users/lucassilva/chromhmm/sra_dir/peaks_hmm/H3K4me1/H3K4me1.bam'},


                  'tss':{'H3K4ME1':'/work/users/lucassilva/chromhmm/sra_dir/peaks_hmm/H3K4me1/H3K4me1.bam',
                              'H3K4ME3':'/work/users/lucassilva/chromhmm/sra_dir/peaks_hmm/H3K4me1/H3K4me1.bam'}



                  }

    bed_target_hash = {'exon':'/work/users/lucassilva/chromhmm/sra_dir/peaks_hmm/AR_24/AR_24.bed',

                       'tss':'/work/users/lucassilva/chromhmm/sra_dir/peaks_hmm/AR_24/AR_24.bed',
                       }




    #get_groups_normalized_count_in_bed_region_given_a_mark(hash_in_groups_bam =bam_groups , bed_target_hash= bed_target_hash)


    #table_df = pd.read_table('/home/lucas/PycharmProjects/Epigenetic_manager/EpiManager/tabela_grupo.txt', delimiter='\t')

    # print table_df[table_df['Count'].map(type)==pd.]

    #table_df['Mark'] = table_df['Mark'].astype(float)

    # print table_df
    #tips = sns.load_dataset("tips")
    #tips.info()
    # table_df['Count'] = table_df['Count'].astype('float')
    # table_df['Group'] = table_df['Group'].astype('category')
    # table_df['Mark'] = table_df['Mark'].astype('category')
    #print table_df.info()
    #sns.boxplot(x="Group", y="Count", hue="Mark", data=table_df, palette="Set2")

    #plt.show()



if __name__ == '__main__':
    sys.exit(main())
