#!/usr/bin/env python3
# load functions -------------------------------------------------------------------------
import sys
all_args = sys.argv[1:]
source_code_at = all_args[0]
sys.path.append(source_code_at)
from Py.libs import *
from Py.func import *

# Define paths ----------------------------------------------------------------------------------------------------------
input_paths = read_json(all_args[1])
save_at = all_args[2]
if not os.path.exists(save_at):
    os.makedirs(save_at, exist_ok = True)
task_name = all_args[3]
out_paths = input_paths # so that more keys-value pairs are appneded to an existing dictionary 
out_paths['M1_cv_runs_five_fold_res_concatenated'] = f'{save_at}/M1_cv_runs_five_fold_res_concatenated.csv'
out_paths['M1_cv_runs_str_res_concatenated'] = f'{save_at}/M1_cv_runs_str_res_concatenated.csv'
out_paths['M1_scenarios_res_concatenated'] = f'{save_at}/M1_scenarios_res_concatenated.csv'
out_paths["M1_kws_old_new"] = f'{save_at}/M1_kws_old_new.csv'
out_paths['M1_hparams_data_all'] = f'{save_at}/M1_hparams_all.csv'
out_paths['M1_hparams_data_all_list'] = f'{save_at}/M1_hparams_all_list.json'

logs_at = f'{save_at}/{task_name}.log'

logging.basicConfig(filename=logs_at, level=logging.DEBUG, filemode='w')
logger = logging.getLogger(__name__)
logging.info(f'Attempt at {datetime.now().strftime("%d/%m/%Y %H:%M:%S")}')

# load data ----------------------------------------------------------------------------------------------------------
check = read_json(input_paths["pred_dirs_paths"])
check_df = pd.DataFrame(check).transpose()
models = set([f'{x.split("#")[0]}#{x.split("#")[1]}' for x in check_df.index.tolist()])

for i in models:
    files = check_df.iloc[check_df.index.str.contains(i),:].loc[:, "pred_at"].to_list()
    file_status = [exists(f'{x}/output.csv') for x in files]
    files_present =  [x for x, y in zip(files, file_status) if y]
    files_absent =  [x for x, y in zip(files, file_status) if not y]
    write_files = False
    if len(files_present) == len(files):
        logging.info(f'all files present. complete result written for {i}')
        write_files = True
    elif len(files_present) != len(files) and len(files_present) > 0:
        logging.info(f'partial files present. {len(files_absent)} files of {len(files)} are absent. partial result written for {i}')
        files_absent_df = pd.Series(files_absent, dtype="string")
        logging.info(f'paths for missign files for {i} are --------------------')
        files_absent_df.to_csv(logs_at, index = False, mode='a')
        write_files = True
    elif len(files_present) == 0:
        logging.info(f'all files missing. nothing written for {i}')
        write_files = False
    else:
        logging.info(f'all files missing. unaccounted error for {i}. please check')
        write_files = False
    
    if write_files:
        # Load files and write them
        out = None
        for path in files_present:
            read_data = pd.read_csv(f'{path}/output.csv')
            if out is None:
                out = read_data
            else:
                out = pd.concat([out, read_data])
                
        out_paths[i] = f'{save_at}/{i}.csv'
        out.to_csv(f'{out_paths[i]}', index = False)

# for hparams
execute = False
if execute:
    all_params_df = pd.DataFrame(check).transpose().reset_index()
    ids = all_params_df["index"].str.split("#", n = 1, expand = True)
    all_params_df["type"] = ids[0]
    all_params_df["run"] = ids[1]
    
    data_out  = None
    data_out_list = {}
    data_out_list_cv = []
    data_out_list_str = []
    data_out_list_sce = []
    for i in range(len(all_params_df)):
        data_json = read_json(f"{all_params_df['model_at'][i]}/best_params.json")
        desired_order_list = ['CNN_f_fl', 'CNN_ks_fl', 'CNN_ap_fl', 
                              'CNN_num_vl', 
                              'CNN_f_vl_0', 'CNN_ks_vl_0', 'CNN_ap_vl_0', 
                              'CNN_f_vl_1', 'CNN_ks_vl_1', 'CNN_ap_vl_1', 
                              'CNN_f_vl_2', 'CNN_ks_vl_2', 'CNN_ap_vl_2', 
                              'CNN_f_vl_3', 'CNN_ks_vl_3', 'CNN_ap_vl_3',
                              'CNN_num_dl', 
                              'CNN_unit_dl_0', 'CNN_drop_rate_dl_0', 
                              'CNN_unit_dl_1', 'CNN_drop_rate_dl_1', 
                              'CNN_unit_dl_2', 'CNN_drop_rate_dl_2', 
                              'CNN_unit_dl_3', 'CNN_drop_rate_dl_3', 
                              'l_rate', 'beta_val_1', 'beta_val_2']
        data_json_reordered = {k: data_json[k] for k in desired_order_list}
        data_json_reordered["run"] = all_params_df['run'][i]
        if(all_params_df["type"][i] == "cv_runs_five_fold"):
            data_out_list_cv.append(data_json_reordered)
        elif (all_params_df["type"][i] == "cv_runs_str"):
            data_out_list_str.append(data_json_reordered)
        elif (all_params_df["type"][i] == "scenarios"):
            data_out_list_sce.append(data_json_reordered)
        data_df = pd.DataFrame(data_json, index=[all_params_df['run'][i]]).transpose()
        if data_out is None:
            data_out = data_df
        else:
            data_out = pd.concat([data_out, data_df], axis=1)
    data_out_list["cv"] = data_out_list_cv
    data_out_list["str"] = data_out_list_str
    data_out_list["sce"] = data_out_list_sce
    
    data_out_reset_index = data_out.reset_index()
    data_out_melt = pd.melt(data_out_reset_index, id_vars=[data_out_reset_index.columns[0]], value_vars = data_out_reset_index.columns[1:len(data_out_reset_index.columns)])
    
    data_out_melt.to_csv(f'{out_paths["M1_hparams_data_all"]}', index = False)
    write_json(data_out_list, f'{out_paths["M1_hparams_data_all_list"]}')
    
    logging.info(f'Hparam data written for all types')

# Finish off the script ------------------------------------------------------------------------------------------------
with open(f'{all_args[4]}', "w") as fp:   
    json.dump(out_paths, fp)
print(f'{task_name} completed successfully. Paths written at {all_args[4]}')