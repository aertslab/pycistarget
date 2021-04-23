import logging
import pandas as pd
import sys
import glob
import pyranges as pr 
import ray
import os
import subprocess
from pybiomart import Dataset
import shutil

from .utils import *

from IPython.display import HTML

def homer_results(homer_dict, name, results='known'):
    if results == 'known':
        file = os.path.join(homer_dict[name].outdir, 'knownResults.html')
    if results == 'denovo':
        file = os.path.join(homer_dict[name].outdir, 'homerResults.html')
    inplace_change(file, 'width="505" height="50"', 'width="1010" height="200"')
    return HTML(file)

def homer_find_motifs_genome(homer_path: str,
                             region_sets: Dict[str, pr.PyRanges],
                             outdir: str,
                             genome: str,
                             size: str = 'given',
                             mask: bool = True,
                             denovo: bool = False,
                             length: str = '8,10,12',
                             n_cpu: int = 1,
                             meme_path: str = None,
                             meme_collection_path: str = None,
                             cistrome_annotation: List[str] = ['Direct_annot', 'Motif_similarity_annot', 'Orthology_annot', 'Motif_similarity_and_Orthology_annot']):
    # Create logger
    level    = logging.INFO
    format   = '%(asctime)s %(name)-12s %(levelname)-8s %(message)s'
    handlers = [logging.StreamHandler(stream=sys.stdout)]
    logging.basicConfig(level = level, format = format, handlers = handlers)
    log = logging.getLogger('pyCistrome')
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
                                length, 
                                meme_path,
                                meme_collection_path,
                                cistrome_annotation) for name in list(bed_paths.keys())])
    ray.shutdown()
    homer_dict={list(bed_paths.keys())[i]: homer_dict[i] for i in range(len(homer_dict))}
    return homer_dict

@ray.remote
def homer_ray(homer_path: str,
              bed_path: str,
              name: str,
              outdir: str,
              genome: str,
              size: str = 'given',
              mask: bool = True,
              denovo: bool = False,
              length: str = '8,10,12',
              meme_path: str = None,
              meme_collection_path: str = None,
              cistrome_annotation: List[str] = ['Direct_annot', 'Motif_similarity_annot', 'Orthology_annot', 'Motif_similarity_and_Orthology_annot']):
    # Create logger
    level    = logging.INFO
    format   = '%(asctime)s %(name)-12s %(levelname)-8s %(message)s'
    handlers = [logging.StreamHandler(stream=sys.stdout)]
    logging.basicConfig(level = level, format = format, handlers = handlers)
    log = logging.getLogger('pyCistrome')
    
    if os.path.exists(outdir):
        shutil.rmtree(outdir)
    os.mkdir(outdir)
    
    log.info('Running '+ name)
    Homer_res = Homer(homer_path, bed_path, name, outdir, genome, size, mask, denovo, length, meme_path, meme_collection_path, cistrome_annotation)
    log.info(name + ' done!')
    return Homer_res

class Homer(): 
    def __init__(self,
                 homer_path: str,
                 bed_path: str,
                 name: str,
                 outdir: str,
                 genome: str,
                 size: str = 'given',
                 mask: bool = True,
                 denovo: bool = False,
                 length: str = '8,10,12',
                 meme_path: str = None,
                 meme_collection_path: str = None,
                 cistrome_annotation: List[str] = ['Direct_annot', 'Motif_similarity_annot', 'Orthology_annot', 'Motif_similarity_and_Orthology_annot']):
        self.homer_path = homer_path
        self.bed_path = bed_path
        self.genome = genome
        self.outdir = outdir
        self.size = size
        self.len = length
        self.mask = mask
        self.denovo = denovo
        self.name = name
        self.meme_path = meme_path
        self.meme_collection_path = meme_collection_path
        self.cistrome_annotation = cistrome_annotation
        self.known_motifs = None
        self.denovo_motifs = None
        self.known_motif_hits = None
        self.denovo_motif_hits = None
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
           
        log.info("Annotating motifs for " + self.name)     
        self.add_motif_annotation_homer()
        log.info("Finding motif hits for " + self.name)   
        self.find_motif_hits(n_cpu=1)
        log.info("Getting cistromes for " + self.name) 
        self.get_cistromes(self.cistrome_annotation)
        
    def load_known(self):
        known = pd.read_csv(os.path.join(self.outdir, 'knownResults.txt'), sep='\t')
        return known
    
    def load_denovo(self):
        denovo = pd.read_html(os.path.join(self.outdir, 'homerResults.html'), header=0)[0].iloc[:,[7,2,3,4,5,6]]
        denovo.iloc[:,0] = [x.split('More')[0] for x in denovo.iloc[:,0]]
        denovo.to_csv(os.path.join(self.outdir, 'homerResults.txt'), sep='\t', index=False)
        return denovo
    
    def add_motif_annotation_homer(self):
        # Create logger
        level    = logging.INFO
        format   = '%(asctime)s %(name)-12s %(levelname)-8s %(message)s'
        handlers = [logging.StreamHandler(stream=sys.stdout)]
        logging.basicConfig(level = level, format = format, handlers = handlers)
        log = logging.getLogger('pyCistrome')
        
        if self.known_motifs is not None:
            if self.known_motifs.shape[0] != 0:
                log.info('Annotating known motifs')
                # Prepare cistarget annotation
                if 'mm' in self.genome:
                    species = 'mus_musculus'
                if 'dm' in self.genome:
                    species = 'drosophila_melanogaster'
                if 'hg' in self.genome:
                    species = 'homo_sapiens'
                ctx_motif_annotation = load_motif_annotations(species)
                motifs = self.known_motifs
                homer_motifs = 'homer__' + motifs['Consensus'] + '_' + [x.split('(')[0] for x in motifs['Motif Name']]

                motifs['MotifID'] = homer_motifs
                homer_motifs = [x for x in homer_motifs if x in ctx_motif_annotation.index.tolist()]
                ctx_motif_annotation = ctx_motif_annotation.loc[list(set(homer_motifs))].reset_index()

                # Prepare homer annotation
                homer_motif_annotation = pd.read_csv(os.path.abspath(os.path.join(os.path.dirname(self.homer_path), '..', 'motifs/extras/table.txt')),
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
                                   else motifs.Direct_annot.tolist()[x] for x in range(motifs.shape[0])]
                motifs.Direct_annot = motifs.Direct_annot.replace('nan, ', '', regex=True)
                motifs.Direct_annot = motifs.Direct_annot.replace(', nan', '', regex=True)
                motifs = motifs.drop(['MotifID', 'Homer_annot'], axis=1)
                self.known_motifs = motifs
        
        if self.denovo_motifs is not None:
            if self.denovo_motifs.shape[0] != 0:
                if self.meme_path is None:
                    log.info('Parameter meme_path is not provided. Skipping annotation of de novo motifs')
                elif self.meme_collection_path is None:
                    log.info('Parameter meme_collection_path is not provided. Skipping annotation of de novo motifs')
                else:
                    # Find closest match for denovo motifs in cistarget database (as meme)
                    log.info('Comparing de novo motifs with given motif collection with tomtom')
                    homer_motif_paths = glob.glob(os.path.join(self.outdir, 'homerResults', 'motif*[0-9$].motif'))
                    homer_motif_paths = [x for x in homer_motif_paths if 'similar' not in x] 
                    tomtom_pd = pd.concat([tomtom(x, self.meme_path, self.meme_collection_path) for x in homer_motif_paths])
                    ctx_motif_annotation = load_motif_annotations(species)
                    homer_motifs = [x for x in tomtom_pd.iloc[:,1].tolist() if x in ctx_motif_annotation.index.tolist()]
                    ctx_motif_annotation = ctx_motif_annotation.loc[list(set(homer_motifs))].reset_index()
                    ctx_motif_annotation = ctx_motif_annotation.rename(columns={'MotifID': 'Best Match/Tomtom'})
                    # Bind cisTarget annotation
                    tomtom_pd = pd.merge(tomtom_pd, ctx_motif_annotation, on='Best Match/Tomtom', how='left')
                    motifs = pd.merge(self.denovo_motifs, tomtom_pd, on='Best Match/Details', how='left')
                    self.denovo_motifs = motifs

    def find_motif_hits(self, n_cpu=1):
        # Create logger
        level    = logging.INFO
        format   = '%(asctime)s %(name)-12s %(levelname)-8s %(message)s'
        handlers = [logging.StreamHandler(stream=sys.stdout)]
        logging.basicConfig(level = level, format = format, handlers = handlers)
        log = logging.getLogger('pyCistrome')
        
        if self.known_motifs is not None:
            if self.known_motifs.shape[0] != 0:
                # Merge all motifs to file
                log.info('Retrieving enriched regions per known motif')
                if os.path.exists(os.path.join(self.outdir, 'knownResults', 'all_motifs.motif')):
                	os.remove(os.path.join(self.outdir, 'knownResults', 'all_motifs.motif'))
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
            if self.denovo_motifs.shape[0] != 0:
                # Merge all motifs to file
                log.info('Retrieving enriched regions per de novo motif')
                if os.path.exists(os.path.join(self.outdir, 'homerResults', 'all_motifs.motif')):
                	os.remove(os.path.join(self.outdir, 'homerResults', 'all_motifs.motif'))
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
                
    def get_cistromes(self, annotation: List[str] = ['Direct_annot', 'Motif_similarity_annot', 'Orthology_annot', 'Motif_similarity_and_Orthology_annot']):
        if self.known_motif_hits is not None:
            tfs = get_TF_list(self.known_motifs)
            cistrome_dict = {tf: get_cistrome_per_TF(self.known_motif_hits,  get_motifs_per_TF(self.known_motifs, tf, motif_column = 'Motif Name', annotation=annotation)) for tf in tfs}
            self.known_cistromes = cistrome_dict 
        if self.denovo_motif_hits is not None:
            tfs = get_TF_list(self.denovo_motifs)
            cistrome_dict = {tf: get_cistrome_per_TF(self.denovo_motif_hits,  get_motifs_per_TF(self.denovo_motifs, tf, motif_column = 'Best Match/Details', annotation=annotation)) for tf in tfs}
            self.denovo_cistromes = cistrome_dict 
            self.denovo_cistromes.keys()
                