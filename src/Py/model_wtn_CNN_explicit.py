from .libs import *
from .func import CNN_net_fixed, CNN_net_flex, dense_post_concat

# Dev:Model 4 ---------------------------------------------------------------------------------------------------------------------
class tuner(kt.HyperModel):
    def __init__(self, marker_n, ec_n, env_n, tune = False, batch_norm = False):
        self.marker_n = marker_n # markers
        self.env_n = env_n # environment variables
        self.ec_n = ec_n # cofactors calculated per environment variable
        self.tune = tune
        self.batch_norm = batch_norm
        print(f'markers = {self.marker_n}, weather_variables = {self.env_n} and ec\'s = {self.ec_n}')
    
    def build(self, hp):
        X_list = [CNN_net_fixed(hp, n_in = self.ec_n, 
                               model_name = f"ec_{x}", 
                               kernel_size = 2) for x in range(self.env_n)]
        g_net_in, g_net_out = CNN_net_flex(hp, n_in = self.marker_n, 
                                        model_name = "g_a",
                                        base_layer = 6, max_layer = 4, tuning = self.tune)
        input_list = [x[0] for x in X_list]
        input_list.append(g_net_in)
        
        X_list_2 = [CNN_net_flex(hp, ec_layer = X_list[x][1], g_layer = g_net_out, model_name = f"pen_{x}", base_layer = 7, max_layer = 3, kernel_size = 2, pool_size = 8, strides = 6) for x in range(len(X_list))]
        
        output_list = [x[1] for x in X_list_2]
        
        model_concat = concatenate(output_list, name = "concat_in")
        
        # do normalization
        model_norm = BatchNormalization()(model_concat)
        
        model_compiled = dense_post_concat(hp, inputs = input_list, concatenated_model = model_norm, base_layer = 5, max_layer = 8)

        return(model_compiled)