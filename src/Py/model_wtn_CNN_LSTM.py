from .libs import *
from .func import LSTM_net, CNN_net_flex, dense_post_concat
# Model 2 ---------------------------------------------------------------------------------------------------------------------
class tuner(kt.HyperModel):
    def __init__(self, marker_n, env_n, days_n, tune = False, batch_norm = False):
        self.marker_n = marker_n # markers
        self.env_n = env_n # environment variables
        self.days_n = days_n # days for which the variable was recorded
        self.tune = tune
        self.batch_norm = batch_norm
        print(f'markers = {self.marker_n}, weather_variables = {self.env_n} and days = {self.days_n}')
        
    def build(self, hp):
        #  LSTM model ------------------------------------------------------------------------------------
        net_env = LSTM_net(hp, n_in_1 = self.days_n, 
                            n_in_2 = self.env_n)
       
        # CNN model -------------------------------------------------------------------------------------- 
        net_g_a = CNN_net_flex(hp, n_in = self.marker_n, 
                          model_name = "g_a", base_layer = 7, max_layer = 4, tuning = self.tune)
        
        # concatenate models -----------------------------------------------------------------------------
        model_concat = concatenate([net_env[1], net_g_a[1]], name="concat_in")
        model_compiled = dense_post_concat(hp, [net_env[0], net_g_a[0]], model_concat)
        
        return model_compiled