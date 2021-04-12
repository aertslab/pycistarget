import logging
import pandas as pd
import sys
import pyranges as pr # with pyfaidx
import ray
import os
import subprocess
from .utils import *

from IPython.display import HTML

def homer_results(homer_dict, name, results='known'):
    if results == 'known':
        file = os.path.join(homer_dict[name].outdir, 'knownResults.html')
    if results == 'denovo':
        file = os.path.join(homer_dict[name].outdir, 'homerResults.html')
    inplace_change(file, 'width="505" height="50"', 'width="1010" height="200"')
    return HTML(file)

def homer_find_motifs_genome(homer_path, region_sets, outdir, genome, size='given', mask=True, denovo=False, length='8,10,12', n_cpu=5):
    # Create logger
    level    = logging.INFO
    format   = '%(asctime)s %(name)-12s %(levelname)-8s %(message)s'
    handlers = [logging.StreamHandler(stream=sys.stdout)]
    logging.basicConfig(level = level, format = format, handlers = handlers)
    log = logging.getLogger('cisTopic')
    # Save regions in dict to the output dir
    bed_paths={}
    bed_dir = os.path.join(outdir, 'regions_bed')
    # Create bed directory
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    if not os.path.exists(bed_dir):
        os.mkdir(bed_dir)
    # Create bed files for Homer
    for key in region_sets.keys():
        bed_path = os.path.join(bed_dir, key+'.bed')
        region_sets[key].Name =  coord_to_region_names(region_sets[key])
        region_sets[key].to_bed(path=bed_path, keep=False, compression='infer', chain=False)
        bed_paths[key] = bed_path
    # Run Homer
    ray.init(num_cpus=n_cpu)
    homer_dict = ray.get([homer_ray.remote(homer_path,
                                bed_paths[name],
                                name,
                                outdir + name, 
                                genome,
                                size,
                                mask,
                                denovo,
                                length) for name in list(bed_paths.keys())])
    ray.shutdown()
    homer_dict={list(bed_paths.keys())[i]: homer_dict[i] for i in range(len(homer_dict))}
    return homer_dict

@ray.remote
def homer_ray(homer_path, bed_path, name, outdir, genome, size='given', mask=True, denovo=False, length='8,10,12'):
    # Create logger
    level    = logging.INFO
    format   = '%(asctime)s %(name)-12s %(levelname)-8s %(message)s'
    handlers = [logging.StreamHandler(stream=sys.stdout)]
    logging.basicConfig(level = level, format = format, handlers = handlers)
    log = logging.getLogger('cisTopic')
    
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    
    log.info('Running '+ name)
    Homer_res = Homer(homer_path, bed_path, name, outdir, genome, size, mask, denovo, length)
    log.info(name + ' done!')
    return Homer_res

class Homer(): 
    def __init__(self, homer_path, bed_path, name, outdir, genome, size='given', mask=True, denovo=False, length='8,10,12'):
        self.homer_path = homer_path
        self.bed_path = bed_path
        self.genome = genome
        self.outdir = outdir
        self.size = size
        self.len = length
        self.mask = mask
        self.denovo = denovo
        self.name = name
        self.known_motifs = None
        self.denovo_motifs = None
        self.known_cistromes = None
        self.denovo_cistromes = None
        self.run()

    def run(self):
        # Create logger
        level    = logging.INFO
        format   = '%(asctime)s %(name)-12s %(levelname)-8s %(message)s'
        handlers = [logging.StreamHandler(stream=sys.stdout)]
        logging.basicConfig(level = level, format = format, handlers = handlers)
        log = logging.getLogger('cisTopic')
        
        if self.mask == True and self.denovo == False:
            cmd = os.path.join(self.homer_path, 'findMotifsGenome.pl') + ' %s %s %s -preparsedDir %s -size %s -len %s -mask -nomotif -keepFiles'
        if self.mask == True and self.denovo == True:
            cmd = os.path.join(self.homer_path, 'findMotifsGenome.pl') + ' %s %s %s -preparsedDir %s -size %s -len %s -mask -keepFiles'
        if self.mask == False and self.denovo == False:
            cmd = os.path.join(self.homer_path, 'findMotifsGenome.pl') + ' %s %s %s -preparsedDir %s -size %s -len %s -nomotif -keepFiles'
        if self.mask == False and self.denovo == True:
            cmd = os.path.join(self.homer_path, 'findMotifsGenome.pl') + ' %s %s %s -preparsedDir %s -size %s -len %s -nomotif -keepFiles'
            
        cmd = cmd % (self.bed_path, self.genome, self.outdir, self.outdir, self.size, self.len)
        log.info("Running Homer for " + self.name + " with %s", cmd)
        try:
            subprocess.check_output(args=cmd, shell=True, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            raise RuntimeError("command '{}' return with error (code {}): {}".format(e.cmd, e.returncode, e.output))
        
        try:
            self.known_motifs = self.load_known()
        except:
            log.info('No known results found')
        if self.denovo == True:
            try:
                self.denovo_motifs = self.load_denovo()
            except:
                log.info('No de novo results found')
        
    def load_known(self):
        known = pd.read_csv(os.path.join(self.outdir, 'knownResults.txt'), sep='\t')
        return known
    
    def load_denovo(self):
        denovo = pd.read_html(os.path.join(self.outdir, 'homerResults.html'), header=0)[0].iloc[:,[7,2,3,4,5,6]]
        denovo.iloc[:,0] = [x.split('More')[0] for x in denovo.iloc[:,0]]
        denovo.to_csv(os.path.join(self.outdir, 'homerResults.txt'), sep='\t', index=False)
        return denovo

    def find_motif_hits(self, n_cpu=1):
        if self.known_motifs is not None:
            # Merge all motifs to file
            for f in glob.glob(os.path.join(self.outdir, 'knownResults', '*.motif')):
                os.system("cat "+f+" >> "+os.path.join(self.outdir, 'knownResults', 'all_motifs.motif'))
            cmd = os.path.join(self.homer_path, 'homer2 find') + ' -s %s -m %s -o %s -p %s'
            cmd = cmd % (os.path.join(self.outdir, 'targetgiven.seq'), os.path.join(self.outdir, 'knownResults', 'all_motifs.motif'), os.path.join(self.outdir, 'knownResults_motif_hits.bed'), n_cpu)
            try:
                subprocess.check_output(args=cmd, shell=True, stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as e:
                raise RuntimeError("command '{}' return with error (code {}): {}".format(e.cmd, e.returncode, e.output))

            known_motif_hits = pd.read_csv(os.path.join(self.outdir, 'knownResults_motif_hits.bed'), sep='\t', header=None)
            self.known_motif_hits = known_motif_hits.groupby(3)[0].apply(lambda g: g.values.tolist()).to_dict()

        if self.denovo_motifs is not None:
            # Merge all motifs to file
            for f in glob.glob(os.path.join(self.outdir, 'homerResults', '*.motif')):
                os.system("cat "+f+" >> "+os.path.join(self.outdir, 'homerResults', 'all_motifs.motif'))
            os.system("sed -i 's/\t.*BestGuess:/\t/g' "+os.path.join(self.outdir, 'homerResults', 'all_motifs.motif'))
            cmd = os.path.join(self.homer_path, 'homer2 find') + ' -s %s -m %s -o %s -p %s'
            cmd = cmd % (os.path.join(self.outdir, 'targetgiven.seq'), os.path.join(self.outdir, 'homerResults', 'all_motifs.motif'), os.path.join(self.outdir, 'homerResults_motif_hits.bed'), n_cpu)
            try:
                subprocess.check_output(args=cmd, shell=True, stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as e:
                raise RuntimeError("command '{}' return with error (code {}): {}".format(e.cmd, e.returncode, e.output))

            denovo_motif_hits = pd.read_csv(os.path.join(self.outdir, 'homerResults_motif_hits.bed'), sep='\t', header=None)
            denovo_motif_hits  = denovo_motif_hits.groupby(3)[0].apply(lambda g: g.values.tolist()).to_dict()
            self.denovo_motif_hits = {k:denovo_motif_hits[k] for k in denovo_motif_hits.keys() if not k[0].isdigit()}
            
    def add_motif_annotation_homer(self,
                              species: str):
        # Prepare cistarget annotation
        ctx_motif_annotation = load_motif_annotations(species)
        homer_motifs = 'homer__' + self.known_motifs['Consensus'] + '_' + [x.split('(')[0] for x in self.known_motifs['Motif Name']]

        motifs['MotifID'] = homer_motifs
        homer_motifs = [x for x in homer_motifs if x in ctx_motif_annotation.index.tolist()]
        ctx_motif_annotation = ctx_motif_annotation.loc[list(set(homer_motifs))].reset_index()

        # Prepare homer annotation
        homer_motif_annotation = pd.read_csv('/data/leuven/software/biomed/haswell_centos7/2018a/software/HOMER/4.10.4-foss-2018a/motifs/extras/table.txt',
                                            sep='\t', error_bad_lines=False).iloc[:,[1,11]].dropna()
        # If not human, convert by homology
        if species is not 'homo_sapiens':
            dataset = Dataset(name='hsapiens_gene_ensembl',
                          host='http://www.ensembl.org')
            if species is 'mus_musculus':
                biomart_query = 'mmusculus_homolog_associated_gene_name'
            if species is 'drosophila_melanogaster':
                biomart_query = 'dmelanogaster_homolog_associated_gene_name'

            human2specie = dataset.query(attributes=['external_gene_name', biomart_query])
            human2specie.index = human2specie['Gene name']
            # Check that the TF has homolog
            TF_names = [x for x in homer_motif_annotation.iloc[:,1].tolist() if x in human2specie.index.tolist()]
            human2specie = human2specie.loc[TF_names,:]
            human2specie.columns = ['Symbol', 'Homolog']
            df = pd.merge(homer_motif_annotation, human2specie, on='Symbol', how='left')
            homer_motif_annotation = df.iloc[:,[0,2]]
            homer_motif_annotation.columns = ['Name', 'Symbol']

        # We first bind the cisTarget annotation
        motifs = pd.merge(motifs, ctx_motif_annotation, on='MotifID', how='left')
        # We now bind the Homer annotation
        homer_motif_annotation.columns = ['Motif Name', 'Homer_annot']
        motifs = pd.merge(motifs, homer_motif_annotation, on='Motif Name', how='left')
        # If Homer_annot is not in Direct_annot we will add it
        # Concatenate
        motifs.Direct_annot = [str(motifs.Direct_annot.tolist()[x]) + ', ' + str(motifs.Homer_annot.tolist()[x]) 
                           if (str(motifs.Homer_annot.tolist()[x]) not in str(motifs.Direct_annot.tolist()[x]))
                           else motifs.Direct_annot.tolist()[x] for x in range(df.shape[0])]
        motifs.Direct_annot = motifs.Direct_annot.replace('nan, ', '', regex=True)
        motifs.Direct_annot = motifs.Direct_annot.replace(', nan', '', regex=True)
        motifs.drop(['MotifID', 'Homer_annot'], axis=1)
        self.known_motifs = motifs