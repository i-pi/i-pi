#import reweight object
from reweight import Reweight
#load the configure file as dict
config_file = "./config.json"
water = Reweight(config_file)
water.run()
