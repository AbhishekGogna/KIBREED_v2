#!/usr/bin/env python3
# load functions -------------------------------------------------------------------------
import sys
source_code_at = '/proj/src'
sys.path.append(source_code_at)
from Py.libs import *
from Py.func import read_json, write_json, write_pkl

proj_paths_raw = read_json("/proj/results/core_paths.json")
proj_paths = {}
for key in proj_paths_raw.keys():
    proj_paths[key] = proj_paths_raw[key][0] # ToDo: this issue is from a josn file written in R. Fix it
#ext_dir_KIBREED = "/proj/ext_dir/KIBREED/results_plots"
res_at = "/proj/results/Py"
res_at_KDG = f"{res_at}/KIBREED_data_generation"
res_at_PCD = f"{res_at}/process_cgm_data"
cv_from_R = "/proj/results/R/generate_prediction_data"

# Define input_files ----------------------------------------------------------------------------------------------------------
input_paths = {}
input_paths['g_data_a'] = f'{res_at_KDG}/GNdata_comb_add.feather'
input_paths['g_data_d'] = f'{res_at_KDG}/GNdata_comb_dom.feather'
#input_paths['connect'] = f'{res_at_KDG}/KIBREED_pheno_data_to_genodata_transition_v1.txt'
input_paths['p_data_BLUEs'] = f'{res_at_KDG}/BLUES_acr_env.feather'
input_paths['p_data_wtn'] = f'{res_at_KDG}/BLUES_within_env.feather'
input_paths['p_data_raw'] = f'{res_at_PCD}/BLUEs_within_env_cgm.feather'
input_paths['e_data_raw'] = f'{res_at_KDG}/climate_data.feather'
input_paths['ec_data'] = f'{res_at_KDG}/ec_mat.feather'
input_paths['cgm_g_s_mat'] = f'{res_at_PCD}/g_s_mat.feather'
input_paths['clust_data'] = f"/proj/results/R/feature_importance/clust_data.csv"

## cv scenarios
input_paths['cv_acr_5f'] = f'{cv_from_R}/cv_acr_5f/cv_acr_5f.json'
input_paths['cv_acr_str'] = f'{cv_from_R}/cv_acr_str/cv_acr_str.json'
input_paths['cv_acr_sce'] = f'{cv_from_R}/cv_acr_sce/cv_acr_sce.json'
input_paths['cv_wtn_tra'] = f'{cv_from_R}/cv_wtn_tra/cv_wtn_tra.json'
input_paths['cv_wtn_cvL'] = f'{cv_from_R}/cv_wtn_cvL/cv_wtn_cvL.json'
input_paths['cv_wtn_LoO'] = f'{cv_from_R}/cv_wtn_LoO/cv_wtn_LoO.json'


## Dump input files
if not os.path.exists(res_at):
    os.makedirs(res_at, exist_ok = True)
write_json(input_paths, f'{res_at}/input_paths.json')

# Define tasks ----------------------------------------------------------------------------------------------------------------
def task_preprocessing():
    '''performs preprocessing'''
    task_name = 'preprocessing'
    input_file_path = f'{res_at}/input_paths.json'
    tag_name = f'{proj_paths["run_Py"]}/{task_name}'
    output_path = f'{proj_paths["results_Py"]}/{task_name}'
    target_file = f'{proj_paths["results_Py"]}/{task_name}_res.json'
    return {
        'file_dep': [input_file_path, f'{tag_name}.py'],
        'targets': [target_file],
        'actions': [f'python3 {tag_name}.py {source_code_at} {input_file_path} {output_path} {task_name} {target_file} > {tag_name}.log 2> {tag_name}.err'],
    }

def task_create_slurm_scripts():
    '''creates slurm scripts'''
    preceding_task = "preprocessing"
    input_file_path = f'{proj_paths["results_Py"]}/{preceding_task}_res.json'
    
    task_name = 'create_slurm_scripts'
    tag_name = f'{proj_paths["run_Py"]}/{task_name}'
    output_path = f'{proj_paths["results_Py"]}/{task_name}'
    target_file = f'{proj_paths["results_Py"]}/{task_name}_res.json'  
    return {
        'file_dep': [input_file_path, f'{tag_name}.py'],
        'targets': [target_file],
        'actions': [f'python3 {tag_name}.py {source_code_at} {input_file_path} {output_path} {task_name} {target_file} > {tag_name}.log 2> {tag_name}.err'],
    }

def task_feature_importance():
    '''detects feature importance'''
    input_file_path = f'{res_at}/input_paths.json'
    
    task_name = 'feature_importance'
    tag_name = f'{proj_paths["run_Py"]}/{task_name}'
    output_path = f'{proj_paths["results_Py"]}/{task_name}'
    target_file = f'{proj_paths["results_Py"]}/{task_name}_res.json'  
    return {
        'file_dep': [input_file_path, f'{tag_name}.py'],
        'targets': [target_file],
        'actions': [f'python3 {tag_name}.py {source_code_at} {input_file_path} {output_path} {task_name} {target_file} > {tag_name}.log 2> {tag_name}.err'],
    }

#def task_submit_jobs():
#    '''submits master scripts to ipk slurm distribution'''
#    preceding_task = "create_slurm_scripts"
#    input_file_path = f'{proj_paths["results_Py"]}/{preceding_task}_res.json'
#
#    task_name = 'submit_jobs'
#    tag_name = f'{proj_paths["run_Py"]}/{task_name}'
#    output_path = f'{proj_paths["results_Py"]}/{task_name}'
#    target_file = f'{proj_paths["results_Py"]}/{task_name}_res.json'  
#    return {
#        'file_dep': [input_file_path, f'{tag_name}.py'],
#        'targets': [target_file],
#        'actions': [f'python3 {tag_name}.py {source_code_at} {input_file_path} {output_path} {task_name} {target_file} > {tag_name}.log 2> {tag_name}.err'],
#    } # figure out a way to do this
#
#def task_process_prediction_output():
#    '''processes the output of the predictions'''
#    preceding_task = "submit_jobs"
#    input_file_path = f'{proj_paths["results_Py"]}/{preceding_task}_res.json'
#    
#    task_name = 'process_prediction_output'
#    tag_name = f'{proj_paths["run_Py"]}/{task_name}'
#    output_path = f'{proj_paths["results_Py"]}/{task_name}'
#    target_file = f'{proj_paths["results_Py"]}/{task_name}_res.json'  
#    return {
#        'file_dep': [input_file_path, f'{tag_name}.py'],
#        'targets': [target_file],
#        'actions': [f'python3 {tag_name}.py {source_code_at} {input_file_path} {output_path} {task_name} {target_file} > {tag_name}.log 2> {tag_name}.err'],
#    }

