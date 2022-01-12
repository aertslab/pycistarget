import re
from attr import attr, has
import h5py
import numpy as np
import pandas as pd
from typing import Mapping, Union
import pyranges as pr

from .motif_enrichment_cistarget import cisTarget
from .motif_enrichment_dem import DEM
from .motif_enrichment_homer import Homer

from copy import deepcopy

from anndata._io.h5ad import write_attribute, read_attribute

TO_SAVE_AS_ATTRIBUTE = [
    'annotation',
    'motif_annotation',
    'specie',
    'annotation_version',
    'path_to_motifs',
    'name',
    'max_bg_regions',
    'adjpval_thr',
    'log2fc_thr',
    'mean_fg_thr',
    'motif_hit_thr',
    'n_cpu',
    'cluster_buster_path',
    'path_to_genome_fasta',
    'path_to_motif_annotations',
    'tmp_dir',
    'homer_path',
    'bed_path',
    'genome',
    'outdir',
    'len',
    'size',
    'mask',
    'denovo',
    'meme_path',
    'meme_collection_path',
    'cistrome_annotation',
    'auc_threshold',
    'nes_threshold',
    'rank_threshold',
    'promoter_space',
    "contrasts",

]


def _writer(
        result: Union[cisTarget, DEM, Homer],
        hdf5_grp_or_file: Union[h5py.File, h5py.Group]):

    if isinstance(hdf5_grp_or_file, h5py.File):
        grp = hdf5_grp_or_file.create_group(name=result.name)
    elif isinstance(hdf5_grp_or_file, h5py.Group):
        grp = hdf5_grp_or_file
    else:
        raise ValueError(
            f'hdf5_grp_or_file should be an instance of h5py.File or h5py.Group, not: {type(hdf5_grp_or_file)}')

    for attribute_name in result.__dir__():
        if attribute_name.startswith('__') and attribute_name.endswith('__'):
            pass
        elif hasattr(getattr(result, attribute_name), '__call__'):
            pass
        elif getattr(result, attribute_name) is None:
            pass
        else:
            if attribute_name in TO_SAVE_AS_ATTRIBUTE or type(getattr(result, attribute_name)) == str or not hasattr(getattr(result, attribute_name), '__iter__'):
                # save strings and non-itterables as attributes
                grp.attrs[attribute_name] = getattr(result, attribute_name)
            else:
                attribute = deepcopy(getattr(result, attribute_name))

                # convert pyranges to pandas.DataFrame
                if isinstance(attribute, pr.PyRanges):
                    attribute = attribute.df

                if isinstance(attribute, pd.DataFrame):
                    if attribute.isna().sum().sum() > 0:
                        attribute = attribute.fillna('nan')

                if hasattr(attribute, '__iter__'):
                    if hasattr(attribute, 'keys'):
                        for key in attribute.keys():
                            if isinstance(attribute[key], pr.PyRanges):
                                attribute[key] = attribute[key].df
                            if isinstance(attribute[key], pd.DataFrame):
                                if attribute[key].isna().sum().sum() > 0:
                                    attribute[key] = attribute[key].fillna(
                                        'nan')
                    else:
                        new_attribute = []
                        for x in attribute:
                            if isinstance(x, pr.PyRanges):
                                new_attribute.append(x.df)
                            if isinstance(x, pd.DataFrame):
                                if x.isna().sum().sum() > 0:
                                    new_attribute.append(x.fillna('nan'))
                        attribute = new_attribute

                write_attribute(grp, attribute_name, attribute, dataset_kwargs = {'compression': 'gzip', 'compression_opts': 9})


def _add_analysis_type_to_grp(grp, result):
    if isinstance(result, cisTarget):
        grp.attrs['analysis_type'] = 'cisTarget'
    elif isinstance(result, DEM):
        grp.attrs['analysis_type'] = 'DEM'
    elif isinstance(result, Homer):
        grp.attrs['analysis_type'] = 'Homer'


def write_hdf5(
        result: Union[Mapping[str, Union[cisTarget, DEM, Homer]], cisTarget, DEM, Homer],
        f_name: str):
    hdf5_file = h5py.File(f_name, 'w')
    try:
        if type(result) == dict:
            if not all([isinstance(x, cisTarget) or isinstance(x, DEM) or isinstance(x, Homer) for x in result.values()]):
                raise ValueError(
                    'result should be a mapping between str keys and instances of either cisTarget, DEM or Homer')
            else:
                for sub_result_name in result.keys():
                    sub_result = result[sub_result_name]
                    grp = hdf5_file.create_group(name=sub_result.name)
                    _writer(sub_result, grp)
                    # add analysis_type attributes, this is important in order to use the correct reader.
                    _add_analysis_type_to_grp(grp, sub_result)

        elif isinstance(result, cisTarget) or isinstance(result, DEM) or isinstance(result, Homer):
            grp = hdf5_file.create_group(name=result.name)
            _writer(result, grp)
            _add_analysis_type_to_grp(grp, result)
        else:
            raise ValueError(
                f'result should be an instance of cisTarget, DEM or Homer or a dict mapping keys to such instances. Not: {type(result)}')
    except Exception as e:
        hdf5_file.close()
        raise(e)
    finally:
        hdf5_file.close()


def _DEM_reader(hdf5_grp: h5py.Group) -> DEM:
    
    hdf5_grp_attributes_dict = dict(hdf5_grp.attrs)

    region_sets = read_attribute( hdf5_grp['region_sets'] )
    region_sets = {
        key: pr.PyRanges(region_sets[key])
        for key in region_sets.keys()}

    DEM_obj = DEM(
        adjpval_thr                     = hdf5_grp_attributes_dict['adjpval_thr']                       if 'adjpval_thr'                    in hdf5_grp_attributes_dict.keys() else 0.05, 
        motif_annotation                = hdf5_grp_attributes_dict['motif_annotation']                  if 'motif_annotation'               in hdf5_grp_attributes_dict.keys() else ['Direct_annot', 'Motif_similarity_annot', 'Orthology_annot', 'Motif_similarity_and_Orthology_annot'], 
        annotation_version              = hdf5_grp_attributes_dict['annotation_version']                if 'annotation_version'             in hdf5_grp_attributes_dict.keys() else 'v9', 
        contrasts                       = hdf5_grp_attributes_dict['contrasts']                         if 'contrasts'                      in hdf5_grp_attributes_dict.keys() else 'Other', 
        log2fc_thr                      = hdf5_grp_attributes_dict['log2fc_thr']                        if 'log2fc_thr'                     in hdf5_grp_attributes_dict.keys() else 1, 
        max_bg_regions                  = hdf5_grp_attributes_dict['max_bg_regions']                    if 'max_bg_regions'                 in hdf5_grp_attributes_dict.keys() else None,
        mean_fg_thr                     = hdf5_grp_attributes_dict['mean_fg_thr']                       if 'mean_fg_thr'                    in hdf5_grp_attributes_dict.keys() else 0, 
        n_cpu                           = hdf5_grp_attributes_dict['n_cpu']                             if 'n_cpu'                          in hdf5_grp_attributes_dict.keys() else 1, 
        name                            = hdf5_grp_attributes_dict['name']                              if 'name'                           in hdf5_grp_attributes_dict.keys() else 'DEM', 
        specie                          = hdf5_grp_attributes_dict['specie']                            if 'specie'                         in hdf5_grp_attributes_dict.keys() else None, 
        tmp_dir                         = hdf5_grp_attributes_dict['tmp_dir']                           if 'tmp_dir'                        in hdf5_grp_attributes_dict.keys() else None,
        dem_db_fname                    = hdf5_grp_attributes_dict['dem_db']                            if 'dem_db'                         in hdf5_grp_attributes_dict.keys() else None,
        motif_hit_thr                   = hdf5_grp_attributes_dict['motif_hit_thr']                     if 'motif_hit_thr'                  in hdf5_grp_attributes_dict.keys() else None,
        fraction_overlap                = hdf5_grp_attributes_dict['fraction_overlap']                  if 'fraction_overlap'               in hdf5_grp_attributes_dict.keys() else 0.4,
        cluster_buster_path             = hdf5_grp_attributes_dict['cluster_buster_path']               if 'cluster_buster_path'            in hdf5_grp_attributes_dict.keys() else None, 
        path_to_genome_fasta            = hdf5_grp_attributes_dict['path_to_genome_fasta']              if 'path_to_genome_fasta'           in hdf5_grp_attributes_dict.keys() else None, 
        path_to_motifs                  = hdf5_grp_attributes_dict['path_to_motifs']                    if 'path_to_motifs'                 in hdf5_grp_attributes_dict.keys() else None, 
        genome_annotation               = hdf5_grp_attributes_dict['genome_annotation']                 if 'genome_annotation'              in hdf5_grp_attributes_dict.keys() else None, 
        promoter_space                  = hdf5_grp_attributes_dict['promoter_space']                    if 'promoter_space'                 in hdf5_grp_attributes_dict.keys() else 1000, 
        path_to_motif_annotations       = hdf5_grp_attributes_dict['path_to_motif_annotations']         if 'path_to_motif_annotations'      in hdf5_grp_attributes_dict.keys() else None, 
        motif_similarity_fdr            = hdf5_grp_attributes_dict['motif_similarity_fdr']              if 'motif_similarity_fdr'           in hdf5_grp_attributes_dict.keys() else 0.001, 
        orthologous_identity_threshold  = hdf5_grp_attributes_dict['orthologous_identity_threshold']    if 'orthologous_identity_threshold' in hdf5_grp_attributes_dict.keys() else 0,
        region_sets = region_sets)

    if 'motif_annotation' in hdf5_grp.keys():
        setattr(DEM_obj, 'motif_annotation', read_attribute(hdf5_grp['cistromes'])) 
    if 'motif_hits' in hdf5_grp.keys():
        setattr(DEM_obj, 'motif_hits', read_attribute(hdf5_grp['motif_hits']))

    if 'motif_enrichment' in hdf5_grp.keys():
        motif_enrichment = read_attribute(hdf5_grp['motif_enrichment'])
        for key in motif_enrichment.keys():
            motif_enrichment[key] = motif_enrichment[key].replace('nan', np.nan)
        setattr(DEM_obj, 'motif_enrichment', motif_enrichment)

    if 'regions_to_db' in hdf5_grp.keys():
        regions_to_db = read_attribute(hdf5_grp['regions_to_db'])
        for key in regions_to_db.keys():
            regions_to_db[key] = regions_to_db[key].replace('nan', np.nan)
        setattr(DEM_obj, 'regions_to_db', regions_to_db)

    return DEM_obj
    

def _Homer_reader(hdf5_grp: h5py.Group) -> Homer:

    hdf5_grp_attributes_dict = dict(hdf5_grp.attrs)

    Homer_obj = Homer(
        homer_path                      = hdf5_grp_attributes_dict['homer_path']                        if 'homer_path'                     in hdf5_grp_attributes_dict.keys() else None, 
        bed_path                        = hdf5_grp_attributes_dict['bed_path']                          if 'bed_path'                       in hdf5_grp_attributes_dict.keys() else None, 
        name                            = hdf5_grp_attributes_dict['name']                              if 'name'                           in hdf5_grp_attributes_dict.keys() else None, 
        outdir                          = hdf5_grp_attributes_dict['outdir']                            if 'outdir'                         in hdf5_grp_attributes_dict.keys() else None, 
        genome                          = hdf5_grp_attributes_dict['genome']                            if 'genome'                         in hdf5_grp_attributes_dict.keys() else None, 
        size                            = hdf5_grp_attributes_dict['size']                              if 'size'                           in hdf5_grp_attributes_dict.keys() else 'given', 
        mask                            = hdf5_grp_attributes_dict['mask']                              if 'mask'                           in hdf5_grp_attributes_dict.keys() else True, 
        denovo                          = hdf5_grp_attributes_dict['denovo']                            if 'denovo'                         in hdf5_grp_attributes_dict.keys() else False, 
        length                          = hdf5_grp_attributes_dict['length']                            if 'length'                         in hdf5_grp_attributes_dict.keys() else '8,10,12', 
        meme_path                       = hdf5_grp_attributes_dict['meme_path']                         if 'meme_path'                      in hdf5_grp_attributes_dict.keys() else None, 
        meme_collection_path            = hdf5_grp_attributes_dict['meme_collection_path']              if 'meme_collection_path'           in hdf5_grp_attributes_dict.keys() else None, 
        cistrome_annotation             = hdf5_grp_attributes_dict['cistrome_annotation']               if 'cistrome_annotation'            in hdf5_grp_attributes_dict.keys() else ['Direct_annot', 'Motif_similarity_annot', 'Orthology_annot', 'Motif_similarity_and_Orthology_annot'], 
        motif_similarity_fdr            = hdf5_grp_attributes_dict['motif_similarity_fdr']              if 'motif_similarity_fdr'           in hdf5_grp_attributes_dict.keys() else 0.001, 
        orthologous_identity_threshold  = hdf5_grp_attributes_dict['orthologous_identity_threshold']    if 'orthologous_identity_threshold' in hdf5_grp_attributes_dict.keys() else 0
    )
    if 'denovo_cistromes' in hdf5_grp.keys():
        denovo_cistromes = read_attribute(hdf5_grp['denovo_cistromes'])
        setattr(Homer_obj, 'denovo_cistromes', denovo_cistromes)
    if 'denovo_motif_hits' in hdf5_grp.keys():
        denovo_motif_hits = read_attribute(hdf5_grp['denovo_motif_hits'])
        setattr(Homer_obj, 'denovo_motif_hits',denovo_motif_hits)
    if 'denovo_motifs' in hdf5_grp.keys():
        denovo_motifs = read_attribute(hdf5_grp['denovo_motifs'])
        denovo_motifs.replace('nan', np.nan, inplace = True)
        setattr(Homer_obj, 'denovo_motifs', denovo_motifs)
    if 'known_cistromes' in hdf5_grp.keys():
        known_cistromes = read_attribute(hdf5_grp['known_cistromes'])
        setattr(Homer_obj, 'known_cistromes', known_cistromes)
    if 'known_motif_hits' in hdf5_grp.keys():
        known_motif_hits = read_attribute(hdf5_grp['known_motif_hits'])
        setattr(Homer_obj, 'known_motif_hits', known_motif_hits)
    if 'known_motifs' in hdf5_grp.keys():
        known_motifs = read_attribute(hdf5_grp['known_motifs'])
        known_motifs.replace('nan', np.nan, inplace = True)
        setattr(Homer_obj, 'known_motifs', known_motifs)
    
    return Homer_obj


def _cisTarget_reader(hdf5_grp: h5py.Group) -> cisTarget:
    
    hdf5_grp_attributes_dict = dict(hdf5_grp.attrs)


    region_set = read_attribute( hdf5_grp['region_set'] )
    region_set = pr.PyRanges(region_set)

    cisTarget_obj = cisTarget(
       region_set                       = region_set, 
       name                             = hdf5_grp_attributes_dict['name']                              if 'name'                               in hdf5_grp_attributes_dict.keys() else None,
       specie                           = hdf5_grp_attributes_dict['specie']                            if 'specie'                             in hdf5_grp_attributes_dict.keys() else None, 
       auc_threshold                    = hdf5_grp_attributes_dict['auc_threshold']                     if 'auc_threshold'                      in hdf5_grp_attributes_dict.keys() else 0.005, 
       nes_threshold                    = hdf5_grp_attributes_dict['nes_threshold']                     if 'nes_threshold'                      in hdf5_grp_attributes_dict.keys() else 3, 
       rank_threshold                   = hdf5_grp_attributes_dict['rank_threshold']                    if 'rank_threshold'                     in hdf5_grp_attributes_dict.keys() else 0.05, 
       path_to_motif_annotations        = hdf5_grp_attributes_dict['path_to_motif_annotations']         if 'path_to_motif_annotations'          in hdf5_grp_attributes_dict.keys() else None, 
       annotation_version               = hdf5_grp_attributes_dict['annotation_version']                if 'annotation_version'                 in hdf5_grp_attributes_dict.keys() else 'v9', 
       annotation                       = hdf5_grp_attributes_dict['annotation']                        if 'annotation'                         in hdf5_grp_attributes_dict.keys() else['Direct_annot', 'Motif_similarity_annot', 'Orthology_annot', 'Motif_similarity_and_Orthology_annot'], 
       motif_similarity_fdr             = hdf5_grp_attributes_dict['motif_similarity_fdr']              if 'motif_similarity_fdr'               in hdf5_grp_attributes_dict.keys() else 0.001, 
       orthologous_identity_threshold   = hdf5_grp_attributes_dict['orthologous_identity_threshold']    if 'orthologous_identity_threshold'     in hdf5_grp_attributes_dict.keys() else 0,
       motifs_to_use                    = hdf5_grp_attributes_dict['motifs_to_use']                     if 'motifs_to_use'                      in hdf5_grp_attributes_dict.keys() else None,
    )

    if 'cistromes' in hdf5_grp.keys():
        cistromes = read_attribute(hdf5_grp['cistromes'])
        for key in cistromes.keys():
            for cistrome in cistromes[key].keys():
                cistromes[key][cistrome] = list(cistromes[key][cistrome])
        setattr(cisTarget_obj, 'cistromes', cistromes)
    if 'motif_enrichment' in hdf5_grp.keys():
        motif_enrichment = read_attribute(hdf5_grp['motif_enrichment'])
        motif_enrichment.replace('nan', np.nan, inplace = True)
        setattr(cisTarget_obj, 'motif_enrichment', motif_enrichment)
    if 'motif_hits' in hdf5_grp.keys():    
        motif_hits = read_attribute(hdf5_grp['motif_hits'])
        for key in motif_hits.keys():
            for motif in motif_hits[key].keys():
                motif_hits[key][motif] = list(motif_hits[key][motif])
        setattr(cisTarget_obj, 'motif_hits', motif_hits)
    if 'regions_to_db' in hdf5_grp.keys():
        regions_to_db = read_attribute(hdf5_grp['regions_to_db'])
        regions_to_db.replace('nan', np.nan)
        setattr(cisTarget_obj, 'regions_to_db', regions_to_db)
    
    return cisTarget_obj


def read_h5ad(f_name) -> Union[Mapping[str, Union[cisTarget, DEM, Homer]], cisTarget, DEM, Homer]:
    hdf5_file = h5py.File(f_name, 'r')
    try:
        if len(hdf5_file.keys()) == 1:
            hdf5_grp = hdf5_file[list(hdf5_file.keys())[0]]
            if hdf5_grp.attrs['analysis_type'] == 'cisTarget':
                return_object = _cisTarget_reader(hdf5_grp)
            elif hdf5_grp.attrs['analysis_type'] == 'DEM':
                return_object = _DEM_reader(hdf5_grp)
            elif hdf5_grp.attrs['analysis_type'] == 'Homer':
                return_object = _Homer_reader(hdf5_grp)
            else:
                raise NotImplementedError(f"A reader for {hdf5_grp.attrs['analysis_type']} does not exists, existing readers: cisTarget, DEM, Homer")
        else:
            return_object = {}
            for key in hdf5_file.keys():
                hdf5_grp = hdf5_file[key]
                if hdf5_grp.attrs['analysis_type'] == 'cisTarget':
                    return_sub_object = _cisTarget_reader(hdf5_grp)
                elif hdf5_grp.attrs['analysis_type'] == 'DEM':
                    return_sub_object = _DEM_reader(hdf5_grp)
                elif hdf5_grp.attrs['analysis_type'] == 'Homer':
                    return_sub_object = _Homer_reader(hdf5_grp)
                else:
                    raise NotImplementedError(f"A reader for {hdf5_grp.attrs['analysis_type']} does not exists, existing readers: cisTarget, DEM, Homer")
                return_object[key] = return_sub_object      
    except Exception as e:
        hdf5_file.close()
        raise(e)
    finally:
        hdf5_file.close()
    return return_object

